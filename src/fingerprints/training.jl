#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2026 Jimmy Dubuisson <jimmy@dubuisson.ch>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#

"""
    Document(seed_indices, seed_weights, label)

A single training instance: a sparse seed vector v_k (vertex indices +
corresponding weights) and its class label. Mirrors v1 DF's frequency-weighted
personalization vector restricted to T'(k).
"""
struct Document
    seed_indices::Vector{Int}
    seed_weights::Vector{Float32}
    label::Int
end

"""
    TrainConfig(; lr=1f-3, epochs=50, batch_size=16, weight_decay=0f0,
                  val_every=1, verbose=true, early_stop_patience=0,
                  checkpoint_metric=:val_acc, restore_best=true,
                  checkpoint_path="", seed=42)

Hyperparameters for `train!`. Set `early_stop_patience=0` to disable early
stopping. `weight_decay=0` uses plain Adam; nonzero values use AdamW.

When `restore_best=true`, the model parameters at the epoch with the best
`checkpoint_metric` (one of `:val_acc`, `:val_loss`) are restored at the end
of training. Requires `val_docs` to be passed to `train!`.

When `checkpoint_path != ""`, the best state is also written to that file
on every improvement — so even a crash mid-training leaves the best-so-far
model on disk. The format is selected by file extension: `.jld2` uses JLD2
(required for classifier heads >2 GB; BSON's wire format has a 2 GB
single-document limit), anything else uses BSON. Load back with
`load_ndf_state(path)` (handles both formats), or directly with
`JLD2.load(path)` / `BSON.load(path)` as appropriate.
"""
Base.@kwdef struct TrainConfig
    lr::Float32 = 1f-3
    epochs::Int = 50
    batch_size::Int = 16
    weight_decay::Float32 = 0f0
    val_every::Int = 1
    verbose::Bool = true
    early_stop_patience::Int = 0
    checkpoint_metric::Symbol = :val_acc
    restore_best::Bool = true
    checkpoint_path::String = ""
    # Baseline "best" to beat — set this when warm-starting from a checkpoint
    # so an early epoch with a worse val doesn't overwrite the saved best.
    initial_best_val_acc::Float32 = -Inf32
    initial_best_val_loss::Float32 = Inf32
    seed::Int = 42
    # Device adapter applied to each per-batch tensor (Φ, mask, y) just before
    # the forward pass. Default identity = CPU. Pass `Flux.gpu` when training
    # on GPU so per-batch arrays move to the model's device. The model and Â
    # must already live on that device.
    device::Function = identity
end

# Atomic-ish state writer. Writes to a temp file then renames so readers
# never see a half-written file. Stored as a NamedTuple wrapping the state
# plus metadata about which epoch / metric value it corresponds to.
# Always serializes CPU-backed arrays so checkpoints are portable across
# CPU↔GPU runs.
#
# Format is chosen by extension: `.jld2` → JLD2 (no 2 GB doc-size limit),
# anything else → BSON (compact, but Int32 document-size header overflows
# for >2 GB classifier heads — e.g. flatten readout with 10K classes × 100K
# vocab).
function _save_state(path::AbstractString, state, epoch::Int,
                     val_acc::Real, val_loss::Real)
    mkpath(dirname(abspath(path)))
    tmp = path * ".tmp"
    cpu_state = Flux.cpu(state)
    if endswith(path, ".jld2")
        JLD2.jldsave(tmp; state=cpu_state, epoch=epoch,
                     val_acc=Float32(val_acc), val_loss=Float32(val_loss))
    else
        BSON.bson(tmp, Dict(:state => cpu_state, :epoch => epoch,
                            :val_acc => Float32(val_acc),
                            :val_loss => Float32(val_loss)))
    end
    mv(tmp, path; force=true)
    return nothing
end

"""
    load_ndf_state(path) -> Dict

Load an NDF checkpoint written by `_save_state`. Auto-dispatches on extension
(`.jld2` → JLD2, else BSON). Returns a `Dict{Symbol,Any}` with keys
`:state`, `:epoch`, `:val_acc`, `:val_loss` regardless of underlying format,
so callers can use `loaded[:state]` uniformly.
"""
function load_ndf_state(path::AbstractString)
    if endswith(path, ".jld2")
        raw = JLD2.load(path)  # String-keyed
        return Dict{Symbol,Any}(Symbol(k) => v for (k, v) in raw)
    else
        return BSON.load(path)
    end
end

# Detect classifier output size. Returns 0 when the model has no classifier.
function _n_classes(model::NDF)
    model.classifier === identity && return 0
    layers = model.classifier isa Chain ? model.classifier.layers : (model.classifier,)
    for layer in Iterators.reverse(layers)
        hasproperty(layer, :weight) && return size(layer.weight, 1)
    end
    return 0
end

# Materialize a batch of documents into the (n × d_in × B) tensor + seed mask
# that NDF's batched forward expects.
function _build_batch(docs::AbstractVector{Document},
                      X::AbstractMatrix{Float32},
                      n::Int)
    B = length(docs)
    d_struct = size(X, 2)
    d_in = 1 + d_struct
    Φ = zeros(Float32, n, d_in, B)
    mask = zeros(Bool, n, B)
    for (b, doc) in enumerate(docs)
        for (i, u) in enumerate(doc.seed_indices)
            Φ[u, 1, b] = doc.seed_weights[i]
            mask[u, b] = true
        end
        if d_struct > 0
            @views Φ[:, 2:end, b] .= X
        end
    end
    return Φ, mask
end

# Evaluate loss + accuracy on a held-out set. Splits into batches to bound
# memory; uses testmode so Dropout is off. `X` and `device` are evaluated in
# the same convention as `train!` — `X` is the CPU-side structural feature
# matrix used by `_build_batch`, and `device` moves the per-batch arrays to
# the model's device.
function _evaluate(model::NDF,
                   Â::AbstractMatrix,
                   docs::AbstractVector{Document},
                   X::AbstractMatrix{Float32},
                   batch_size::Int,
                   n_classes::Int;
                   device::Function = identity,
                   central_nodes::Union{Nothing,AbstractVector{<:Integer}}=nothing)
    n = size(X, 1)
    Flux.testmode!(model)
    total_loss = 0.0f0
    correct = 0
    total = 0
    for start in 1:batch_size:length(docs)
        stop = min(start + batch_size - 1, length(docs))
        batch = docs[start:stop]
        Φ_cpu, mask_cpu = _build_batch(batch, X, n)
        labels = Int[d.label for d in batch]
        y_cpu = Flux.onehotbatch(labels, 1:n_classes)
        Φ = device(Φ_cpu); mask = device(mask_cpu); y = device(y_cpu)
        ŷ = model(Φ, Â; seed_mask=mask, central_nodes=central_nodes)
        total_loss += Float32(Flux.logitcrossentropy(ŷ, y)) * length(batch)
        # argmax on possibly-GPU array → move to CPU for comparison with CPU labels.
        ŷ_cpu = Array(ŷ)
        preds = vec(map(I -> I[1], argmax(ŷ_cpu; dims=1)))
        correct += count(preds .== labels)
        total += length(batch)
    end
    Flux.trainmode!(model)
    return total_loss / total, correct / total
end

"""
    train!(model, Â, docs; X=nothing, val_docs=nothing, config=TrainConfig())

Train an `NDF` classifier by minibatch SGD on a list of `Document`s sharing
the same domain graph (represented by its normalized adjacency `Â`).

# Arguments
- `model::NDF`: the model to train in place.
- `Â::AbstractMatrix`: symmetric-normalized adjacency from `prepare_adjacency`.
- `docs::AbstractVector{Document}`: training documents.

# Keyword arguments
- `X`: per-vertex structural feature matrix `(n × d_struct)`. If `nothing`,
       the seed vector alone is used (model must be built with `d_in=1`).
       Otherwise the model must be built with `d_in = 1 + size(X, 2)`.
- `val_docs`: optional held-out documents for periodic evaluation and early
       stopping.
- `config`: `TrainConfig` with hyperparameters.

Returns a `NamedTuple` `(train_loss, train_acc, val_loss, val_acc, test_loss,
test_acc)` where each field is a `Vector{Float32}` of per-evaluation values.
`val_*` are empty when no `val_docs` are supplied. `test_*` are empty when no
`test_docs` are supplied; otherwise they have the same length as `val_*`,
with `NaN` at val evaluations where the checkpoint metric did NOT improve,
and the actual test value at evaluations where it did. This is a pure
diagnostic — `test_docs` is NEVER used for early stopping or checkpoint
selection.
"""
function train!(model::NDF,
                Â::AbstractMatrix,
                docs::AbstractVector{Document};
                X::Union{AbstractMatrix,Nothing}=nothing,
                val_docs::Union{AbstractVector{Document},Nothing}=nothing,
                test_docs::Union{AbstractVector{Document},Nothing}=nothing,
                config::TrainConfig=TrainConfig(),
                central_nodes::Union{Nothing,AbstractVector{<:Integer}}=nothing)
    n = size(Â, 1)
    # X_eff is kept CPU-side; _build_batch reads it to fill structural feature
    # columns before the result is moved to the model's device via `config.device`.
    X_eff::Matrix{Float32} = X === nothing ? zeros(Float32, n, 0) :
                             (X isa Matrix{Float32} ? X : Matrix{Float32}(X))
    dev = config.device

    n_classes = _n_classes(model)
    n_classes > 0 || throw(ArgumentError(
        "train! requires a classifier; build NDF with n_classes > 0"))

    optimizer = config.weight_decay > 0f0 ?
        Flux.AdamW(config.lr, (0.9f0, 0.999f0), config.weight_decay) :
        Flux.Adam(config.lr)
    opt_state = Flux.setup(optimizer, model)

    rng = MersenneTwister(config.seed)
    n_train = length(docs)
    history = (
        train_loss = Float32[],
        train_acc  = Float32[],
        val_loss   = Float32[],
        val_acc    = Float32[],
        test_loss  = Float32[],
        test_acc   = Float32[],
    )
    config.checkpoint_metric ∈ (:val_acc, :val_loss) ||
        throw(ArgumentError("checkpoint_metric must be :val_acc or :val_loss"))

    best_val_loss = config.initial_best_val_loss
    best_val_acc  = config.initial_best_val_acc
    best_state    = nothing
    best_epoch    = 0
    patience      = 0

    for epoch in 1:config.epochs
        Flux.trainmode!(model)
        perm = Random.shuffle(rng, collect(1:n_train))
        epoch_loss = 0.0f0
        epoch_correct = 0

        for start in 1:config.batch_size:n_train
            stop = min(start + config.batch_size - 1, n_train)
            batch_docs = docs[perm[start:stop]]
            Φ_cpu, mask_cpu = _build_batch(batch_docs, X_eff, n)
            labels = Int[d.label for d in batch_docs]
            y_cpu = Flux.onehotbatch(labels, 1:n_classes)
            Φ = dev(Φ_cpu); mask = dev(mask_cpu); y = dev(y_cpu)

            loss, grads = Flux.withgradient(model) do m
                Flux.logitcrossentropy(m(Φ, Â; seed_mask=mask, central_nodes=central_nodes), y)
            end
            Flux.update!(opt_state, model, grads[1])

            B = length(batch_docs)
            epoch_loss += Float32(loss) * B
            # Cheap train accuracy from the same forward (recomputed in
            # testmode would double the cost; we accept the dropout-induced
            # noise here). Move logits to CPU for argmax comparison.
            ŷ_cpu = Array(model(Φ, Â; seed_mask=mask, central_nodes=central_nodes))
            preds = vec(map(I -> I[1], argmax(ŷ_cpu; dims=1)))
            epoch_correct += count(preds .== labels)
        end

        train_loss = epoch_loss / n_train
        train_acc  = epoch_correct / n_train
        push!(history.train_loss, train_loss)
        push!(history.train_acc,  train_acc)

        do_val = val_docs !== nothing && epoch % config.val_every == 0
        if do_val
            val_loss, val_acc = _evaluate(model, Â, val_docs, X_eff,
                                          config.batch_size, n_classes;
                                          device=dev,
                                          central_nodes=central_nodes)
            push!(history.val_loss, val_loss)
            push!(history.val_acc,  val_acc)
            config.verbose && @info "epoch=$epoch train=$(round(train_loss; digits=4))/$(round(train_acc; digits=3)) val=$(round(val_loss; digits=4))/$(round(val_acc; digits=3))"

            # Track the best epoch by the selected metric.
            improved = config.checkpoint_metric === :val_acc ?
                       (val_acc > best_val_acc + 1f-5) :
                       (val_loss < best_val_loss - 1f-5)
            if improved
                best_val_acc  = max(best_val_acc, val_acc)
                best_val_loss = min(best_val_loss, val_loss)
                best_epoch    = epoch
                if config.restore_best
                    best_state = Flux.state(model)
                end
                if !isempty(config.checkpoint_path)
                    _save_state(config.checkpoint_path, best_state === nothing ? Flux.state(model) : best_state, epoch, val_acc, val_loss)
                end
                # Diagnostic test eval at each best epoch — never feeds back
                # into early stop or checkpoint selection.
                if test_docs !== nothing
                    test_loss, test_acc = _evaluate(model, Â, test_docs, X_eff,
                                                    config.batch_size, n_classes;
                                                    device=dev,
                                                    central_nodes=central_nodes)
                    push!(history.test_loss, test_loss)
                    push!(history.test_acc,  test_acc)
                    config.verbose && @info "  test (diagnostic) epoch=$epoch test_loss=$(round(test_loss; digits=4)) test_acc=$(round(test_acc; digits=4))"
                end
                patience = 0
            else
                if test_docs !== nothing
                    push!(history.test_loss, NaN32)
                    push!(history.test_acc,  NaN32)
                end
                patience += 1
                if config.early_stop_patience > 0 &&
                   patience >= config.early_stop_patience
                    config.verbose && @info "early stop at epoch=$epoch (best epoch=$best_epoch val_acc=$(round(best_val_acc; digits=3)) val_loss=$(round(best_val_loss; digits=4)))"
                    break
                end
            end
        elseif config.verbose
            @info "epoch=$epoch train=$(round(train_loss; digits=4))/$(round(train_acc; digits=3))"
        end
    end

    if config.restore_best && best_state !== nothing
        Flux.loadmodel!(model, best_state)
        config.verbose && @info "restored best checkpoint" best_epoch best_val_acc=round(best_val_acc; digits=3) best_val_loss=round(best_val_loss; digits=4)
    end

    return history
end

"""
    fingerprints(model, Â, docs; X=nothing, batch_size=16)

Compute the pooled fingerprint vectors for a collection of documents. Returns
a `(hidden × N)` matrix where each column is one document's fingerprint.
Useful for downstream retrieval/clustering or for plugging the NDF output
into a non-Flux classifier (e.g. AdaBoost+DT, matching v1 DF's evaluation).
"""
function fingerprints(model::NDF,
                      Â::AbstractMatrix,
                      docs::AbstractVector{Document};
                      X::Union{AbstractMatrix,Nothing}=nothing,
                      batch_size::Int=16)
    n = size(Â, 1)
    X_eff::Matrix{Float32} = X === nothing ? zeros(Float32, n, 0) :
                             (X isa Matrix{Float32} ? X : Matrix{Float32}(X))

    # Run model with classifier replaced by identity so we get the pooled
    # fingerprint regardless of whether the trained model has a classifier.
    fp_model = NDF(model.encoder, identity, model.K, model.α, model.readout,
                   model.propagation, model.gamma, model.Wprop)
    Flux.testmode!(fp_model)

    N = length(docs)
    H = nothing
    pos = 1
    for start in 1:batch_size:N
        stop = min(start + batch_size - 1, N)
        Φ, mask = _build_batch(docs[start:stop], X_eff, n)
        f = fp_model(Φ, Â; seed_mask=mask)  # (hidden × B)
        if H === nothing
            H = zeros(Float32, size(f, 1), N)
        end
        H[:, pos:pos + size(f, 2) - 1] .= f
        pos += size(f, 2)
    end
    return H
end
