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

# =============================================================================
# Per-document graph-of-words (TextING-lite).
#
# Motivation. The NDF investigation established that word-only diffusion over a
# FIXED GLOBAL word graph is a reparametrization of bag-of-words and cannot beat
# a BoW-MLP (the multiplicative gate ties WideMLP at Ohsumed 66.1). The signal
# BoW cannot reconstruct is WITHIN-document structure: which words co-occur, in
# which local arrangement, inside THIS document. A per-document graph-of-words
# encodes exactly that — each document becomes its OWN small graph whose nodes
# are the document's words and whose edges are within-window co-occurrences —
# and this family (TextING / TextSSL / TensorGCN) is what actually exceeds BoW
# (~68-70 on Ohsumed).
#
# This file provides the load-bearing infrastructure for that path AND for the
# eventual GLOBAL+PER-DOC fusion (use the corpus graph to denoise / enrich the
# per-doc operator). It is a genuinely different data flow from the global-vocab
# NDF forward: each minibatch is assembled as a BLOCK-DIAGONAL graph over the
# union of its documents' word-nodes, so a single sparse SpMM propagates every
# document in the batch independently. Node features come from a learnable word
# embedding table indexed by GLOBAL vocab id (so embeddings are shared across
# documents — inductive), and a segment (mean + attention) readout pools each
# block back to one document vector.
#
# Reuses `_spmm` (defined in ndf.jl) for every sparse multiply, so the same
# GPU-safe rrule (NoTangent for the operator, cotangent materialized on the
# primal's array backend) covers propagation AND pooling on CUDA.
# =============================================================================

"""
    LocalGraph

A single document rendered as its own graph-of-words.

Fields:
- `node_ids` — GLOBAL vocab id (1..V) of each local node; row `i` of the batch
  embedding gather is `emb[node_ids[i], :]`.
- `src`, `dst` — undirected edge endpoints as LOCAL indices (1..length(node_ids)).
- `ew` — edge weight (within-window co-occurrence count), parallel to `src`/`dst`.
- `seed_w` — per-node seed weight (normalized term frequency), a cheap lexical
  feature concatenated to the embedding.
- `label` — class id.
"""
struct LocalGraph
    node_ids::Vector{Int32}
    src::Vector{Int32}
    dst::Vector{Int32}
    ew::Vector{Float32}
    seed_w::Vector{Float32}
    label::Int
end

"""
    build_local_window_graph(tokens, vocab; window=3, weight=:freq) -> LocalGraph

Build one document's graph-of-words. Nodes are the document's distinct in-vocab
words; an undirected edge connects two words for every time they co-occur within
`window` positions in the token sequence (OOV tokens are dropped, closing the
gap so the window slides over surviving tokens). Edge weights accumulate the
co-occurrence count. `label` is filled by the caller.

`weight` sets the per-node `seed_w` lexical feature: `:freq` = count/|doc|,
`:binary` = 1 for present words. Returns `nothing` when the document has no
in-vocab words (caller drops it).
"""
function build_local_window_graph(tokens::AbstractVector{<:AbstractString},
                                  vocab::Dict{String,Int};
                                  window::Int=3, weight::Symbol=:freq)
    window >= 2 || throw(ArgumentError("window must be >= 2"))
    # Map tokens → global ids, dropping OOV, preserving order.
    seq = Int[]
    for t in tokens
        id = get(vocab, String(t), 0)
        id == 0 && continue
        push!(seq, id)
    end
    isempty(seq) && return nothing

    # Distinct global ids → local node index (first-seen order, deterministic).
    local_of = Dict{Int,Int}()
    node_ids = Int32[]
    counts   = Float32[]
    for g in seq
        li = get(local_of, g, 0)
        if li == 0
            push!(node_ids, Int32(g))
            push!(counts, 1f0)
            local_of[g] = length(node_ids)
        else
            counts[li] += 1f0
        end
    end
    L = length(node_ids)

    # Accumulate within-window co-occurrence counts between DISTINCT local nodes.
    ew = Dict{Tuple{Int32,Int32},Float32}()
    n = length(seq)
    @inbounds for p in 1:n
        gp = local_of[seq[p]]
        qmax = min(n, p + window - 1)
        for q in (p + 1):qmax
            gq = local_of[seq[q]]
            gp == gq && continue
            a = Int32(min(gp, gq)); b = Int32(max(gp, gq))
            ew[(a, b)] = get(ew, (a, b), 0f0) + 1f0
        end
    end

    src = Int32[]; dst = Int32[]; w = Float32[]
    for ((a, b), c) in ew
        push!(src, a); push!(dst, b); push!(w, c)
    end

    total = sum(counts)
    seed_w = weight === :binary ? ones(Float32, L) : (counts ./ total)
    return LocalGraph(node_ids, src, dst, w, seed_w, 0)
end

"""
    docs_to_local_graphs(docs, vocab, label2id; window=3, weight=:freq)
        -> (train::Vector{LocalGraph}, test::Vector{LocalGraph}, n_classes)

Convert a `Vector{TextGcnDoc}` into per-document graphs-of-words, split by the
benchmark's fixed train/test tag. Mirrors `textgcn_to_documents` but emits
`LocalGraph`s instead of global seed vectors. Documents with no in-vocab words
are dropped.
"""
function docs_to_local_graphs(docs::Vector{TextGcnDoc},
                              vocab::Dict{String,Int},
                              label2id::Dict{String,Int};
                              window::Int=3, weight::Symbol=:freq)
    train = LocalGraph[]; test = LocalGraph[]
    for d in docs
        lg = build_local_window_graph(d.tokens, vocab; window=window, weight=weight)
        lg === nothing && continue
        lab = label2id[d.label]
        lg = LocalGraph(lg.node_ids, lg.src, lg.dst, lg.ew, lg.seed_w, lab)
        (d.split === :train ? train : test) |> v -> push!(v, lg)
    end
    return train, test, length(label2id)
end

"""
    LocalBatch

A minibatch of `LocalGraph`s assembled as one block-diagonal graph. All arrays
are CPU-side; the harness moves them to the device (sparse operators become
`CuSparseMatrixCSR`, dense/index arrays become `CuArray`).

Fields:
- `node_ids` — `(N,)` global vocab ids for every node in the batch (for the
  embedding gather), N = Σ block sizes.
- `Â` — `(N × N)` symmetric-normalized block-diagonal adjacency (`prepare_adjacency`).
- `seed_w` — `(N,)` per-node lexical feature.
- `Ind` — `(N × B)` 0/1 node→document indicator (broadcast segment scalars back).
- `Indt` — `(B × N)` its transpose (segment reduction: `Indt * H` sums per doc).
- `sizes` — `(B,)` nodes per document (mean-pool normalizer).
- `labels` — `(B,)` class ids.
"""
# Fields are fully parametric so a batch can be moved to the GPU (node_ids,
# seed_w, sizes → CuArray; Â, Ind, Indt → CuSparseMatrixCSR) by reconstructing
# the same struct. `labels` stays CPU-side (only used for onehot + comparison).
struct LocalBatch{NI,S,SW,SM,SZ}
    node_ids::NI
    Â::S
    seed_w::SW
    Ind::SM
    Indt::SM
    sizes::SZ
    labels::Vector{Int}
end

"""
    build_local_batch(graphs) -> LocalBatch

Assemble a vector of `LocalGraph`s into a block-diagonal `LocalBatch`. Edge
lists are concatenated with per-block offsets so the resulting adjacency is
block-diagonal (documents never share nodes or edges), then symmetric-normalized
with self-loops via `prepare_adjacency`.
"""
function build_local_batch(graphs::AbstractVector{LocalGraph})
    B = length(graphs)
    sizes_i = Int[length(g.node_ids) for g in graphs]
    offsets = cumsum(vcat(0, sizes_i[1:end-1]))
    N = sum(sizes_i)

    node_ids = Vector{Int32}(undef, N)
    seed_w   = Vector{Float32}(undef, N)
    # Edge COO (block-offset, symmetric: push both directions).
    esrc = Int[]; edst = Int[]; ew = Float32[]
    # Segment indicator COO.
    ind_r = Int[]; ind_c = Int[]
    labels = Vector{Int}(undef, B)

    for (b, g) in enumerate(graphs)
        off = offsets[b]
        L = sizes_i[b]
        @inbounds for i in 1:L
            node_ids[off + i] = g.node_ids[i]
            seed_w[off + i]   = g.seed_w[i]
            push!(ind_r, off + i); push!(ind_c, b)
        end
        @inbounds for k in 1:length(g.src)
            a = off + Int(g.src[k]); c = off + Int(g.dst[k])
            push!(esrc, a); push!(edst, c); push!(ew, g.ew[k])
            push!(esrc, c); push!(edst, a); push!(ew, g.ew[k])
        end
        labels[b] = g.label
    end

    A = sparse(esrc, edst, ew, N, N)
    Â = prepare_adjacency(A; self_loops=true)     # D^{-1/2}(A+I)D^{-1/2}, block-diagonal
    Ind  = sparse(ind_r, ind_c, ones(Float32, N), N, B)
    Indt = sparse(ind_c, ind_r, ones(Float32, N), B, N)
    sizes = Float32.(sizes_i)
    return LocalBatch(node_ids, Â, seed_w, Ind, Indt, sizes, labels)
end

# -----------------------------------------------------------------------------
# Model
# -----------------------------------------------------------------------------

"""
    PerDocGNN(V, n_classes; d_emb=200, hidden=200, n_layers=2, dropout=0.5f0,
              readout=:meanattn, head_hidden=256)

TextING-lite graph-of-words classifier. A learnable word-embedding table
(`V × d_emb`, shared across documents = inductive) feeds a stack of `n_layers`
GCN layers run over the block-diagonal batch adjacency, then a segment readout
pools each document's nodes into one vector for a WideMLP-style head.

`readout`:
- `:mean`     — masked mean over the document's nodes.
- `:meanattn` — `[mean ‖ attention-weighted sum]` (segment-softmax attention;
  the TextING readout). Default; concat dim = `2*hidden`.

All matmuls keep the node axis first (`N × d`) and go through `_spmm`, so the
same GPU-safe sparse rrule covers propagation and pooling.
"""
struct PerDocGNN{EM,IP,LW,AT,DR,CL}
    emb::EM        # d_emb × V  (column = word embedding; gathered by node id)
    Win::IP        # (d_emb + 1) × hidden  input projection (embedding ‖ seed_w)
    Ws::LW         # NTuple{n_layers} of hidden × hidden GCN transforms
    att::AT        # hidden × 1  attention scorer (:meanattn only)
    drop::DR       # Dropout (reused across layers)
    classifier::CL # Chain on (readout_dim × B)
    K::Int
    α::Float32
    readout::Symbol
end
Flux.@layer PerDocGNN trainable=(emb, Win, Ws, att, classifier)

function PerDocGNN(V::Int, n_classes::Int;
                   d_emb::Int=200, hidden::Int=200, n_layers::Int=2,
                   dropout::Float32=0.5f0, readout::Symbol=:meanattn,
                   head_hidden::Int=256)
    readout ∈ (:mean, :meanattn) ||
        throw(ArgumentError("readout must be :mean or :meanattn"))
    n_layers >= 1 || throw(ArgumentError("n_layers must be >= 1"))
    n_classes > 0 || throw(ArgumentError("PerDocGNN requires n_classes > 0"))

    glorot(a, b) = Float32.((rand(Float32, a, b) .- 0.5f0) .* (2f0 * sqrt(6f0 / (a + b))))
    emb = glorot(d_emb, V)                 # column per word (NNlib.gather by id)
    Win = glorot(d_emb + 1, hidden)
    Ws  = ntuple(_ -> begin
        W = zeros(Float32, hidden, hidden)
        @inbounds for i in 1:hidden; W[i, i] = 1f0; end   # init ≈ identity propagation
        W
    end, n_layers)
    att = glorot(hidden, 1)
    drop = Flux.Dropout(dropout)
    readout_dim = readout === :meanattn ? 2 * hidden : hidden
    classifier = Chain(Dense(readout_dim => head_hidden, relu),
                       Flux.Dropout(dropout),
                       Dense(head_hidden => n_classes))
    return PerDocGNN(emb, Win, Ws, att, drop, classifier, n_layers, 0.15f0, readout)
end

# Segment-softmax attention pool. `H` is (N × hidden); returns (B × hidden).
# Numerically stabilized by a global scalar shift (segment-invariant).
function _seg_attn_pool(H::AbstractMatrix, att::AbstractMatrix,
                        Ind::AbstractMatrix, Indt::AbstractMatrix)
    s = H * att                                   # (N × 1)
    s = s .- maximum(s)                           # scalar shift (softmax-invariant)
    es = exp.(s)                                  # (N × 1)
    seg = _spmm(Indt, es)                         # (B × 1) per-doc exp sums
    denom = _spmm(Ind, seg)                       # (N × 1) broadcast back to nodes
    w = es ./ denom                               # (N × 1) segment-softmax weights
    Hw = H .* w                                   # (N × hidden)
    return _spmm(Indt, Hw)                        # (B × hidden)
end

"""
    (m::PerDocGNN)(batch::LocalBatch)

Forward pass over one block-diagonal minibatch. Returns logits `(n_classes × B)`.
"""
function (m::PerDocGNN)(node_ids, Â, seed_w, Ind, Indt, sizes)
    x = permutedims(Flux.NNlib.gather(m.emb, node_ids))  # (N × d_emb) GPU-safe gather
    feats = hcat(x, reshape(seed_w, :, 1))        # (N × d_emb+1)
    H = feats * m.Win                             # (N × hidden)
    for W in m.Ws
        H = m.drop(relu.(_spmm(Â, H * W)))        # GCN layer: ReLU(Â (H W))
    end
    mean_pool = _spmm(Indt, H) ./ reshape(sizes, :, 1)   # (B × hidden)
    doc = if m.readout === :meanattn
        hcat(mean_pool, _seg_attn_pool(H, m.att, Ind, Indt))  # (B × 2hidden)
    else
        mean_pool
    end
    return m.classifier(permutedims(doc))         # (n_classes × B)
end

(m::PerDocGNN)(b::LocalBatch) =
    m(b.node_ids, b.Â, b.seed_w, b.Ind, b.Indt, b.sizes)

# -----------------------------------------------------------------------------
# Training
# -----------------------------------------------------------------------------

# Evaluate loss + accuracy over a set of local graphs. `movebatch` moves a
# CPU LocalBatch to the model's device (identity on CPU).
function _evaluate_perdoc(model::PerDocGNN, graphs::AbstractVector{LocalGraph},
                         batch_size::Int, n_classes::Int; movebatch::Function=identity)
    Flux.testmode!(model)
    total_loss = 0f0; correct = 0; total = 0
    for start in 1:batch_size:length(graphs)
        stop = min(start + batch_size - 1, length(graphs))
        cb = build_local_batch(graphs[start:stop])
        b = movebatch(cb)
        y = Flux.onehotbatch(cb.labels, 1:n_classes)
        ŷ = model(b)
        ŷ_cpu = Array(ŷ)
        total_loss += Float32(Flux.logitcrossentropy(ŷ_cpu, y)) * length(cb.labels)
        preds = vec(map(I -> I[1], argmax(ŷ_cpu; dims=1)))
        correct += count(preds .== cb.labels)
        total += length(cb.labels)
    end
    Flux.trainmode!(model)
    return total_loss / total, correct / total
end

"""
    train_perdoc!(model, train_graphs; val_graphs=nothing, test_graphs=nothing,
                  n_classes, config=TrainConfig(), movebatch=identity) -> history

Minibatch SGD for `PerDocGNN`. Mirrors `train!` (AdamW/Adam, per-epoch shuffle,
early stopping + best-checkpoint restore on val accuracy, diagnostic test eval
at each best epoch). `movebatch` moves a CPU `LocalBatch` to the model's device;
pass a CUDA mover for GPU. Returns the same `history` NamedTuple as `train!`.
"""
function train_perdoc!(model::PerDocGNN,
                       train_graphs::AbstractVector{LocalGraph};
                       val_graphs::Union{AbstractVector{LocalGraph},Nothing}=nothing,
                       test_graphs::Union{AbstractVector{LocalGraph},Nothing}=nothing,
                       n_classes::Int,
                       config::TrainConfig=TrainConfig(),
                       movebatch::Function=identity,
                       movearr::Function=identity)
    optimizer = config.weight_decay > 0f0 ?
        Flux.AdamW(config.lr, (0.9f0, 0.999f0), config.weight_decay) :
        Flux.Adam(config.lr)
    opt_state = Flux.setup(optimizer, model)
    rng = MersenneTwister(config.seed)
    n_train = length(train_graphs)
    history = (train_loss=Float32[], train_acc=Float32[], val_loss=Float32[],
               val_acc=Float32[], test_loss=Float32[], test_acc=Float32[])
    best_val_acc = config.initial_best_val_acc
    best_state = nothing; best_epoch = 0; patience = 0

    for epoch in 1:config.epochs
        Flux.trainmode!(model)
        perm = Random.shuffle(rng, collect(1:n_train))
        epoch_loss = 0f0; epoch_correct = 0
        for start in 1:config.batch_size:n_train
            stop = min(start + config.batch_size - 1, n_train)
            batch_graphs = train_graphs[perm[start:stop]]
            cb = build_local_batch(batch_graphs)
            b = movebatch(cb)
            y = movearr(Flux.onehotbatch(cb.labels, 1:n_classes))
            loss, grads = Flux.withgradient(model) do m
                Flux.logitcrossentropy(m(b), y)
            end
            Flux.update!(opt_state, model, grads[1])
            epoch_loss += Float32(loss) * length(cb.labels)
            ŷ_cpu = Array(model(b))
            preds = vec(map(I -> I[1], argmax(ŷ_cpu; dims=1)))
            epoch_correct += count(preds .== cb.labels)
        end
        push!(history.train_loss, epoch_loss / n_train)
        push!(history.train_acc, epoch_correct / n_train)

        if val_graphs !== nothing && epoch % config.val_every == 0
            val_loss, val_acc = _evaluate_perdoc(model, val_graphs,
                config.batch_size, n_classes; movebatch=movebatch)
            push!(history.val_loss, val_loss); push!(history.val_acc, val_acc)
            config.verbose && @info "epoch=$epoch train=$(round(epoch_loss/n_train; digits=4))/$(round(epoch_correct/n_train; digits=3)) val=$(round(val_loss; digits=4))/$(round(val_acc; digits=3))"
            if val_acc > best_val_acc + 1f-5
                best_val_acc = val_acc; best_epoch = epoch
                config.restore_best && (best_state = Flux.state(model))
                if test_graphs !== nothing
                    tl, ta = _evaluate_perdoc(model, test_graphs,
                        config.batch_size, n_classes; movebatch=movebatch)
                    push!(history.test_loss, tl); push!(history.test_acc, ta)
                    config.verbose && @info "  test (diagnostic) epoch=$epoch test_acc=$(round(ta; digits=4))"
                end
                patience = 0
            else
                test_graphs !== nothing && (push!(history.test_loss, NaN32); push!(history.test_acc, NaN32))
                patience += 1
                if config.early_stop_patience > 0 && patience >= config.early_stop_patience
                    config.verbose && @info "early stop at epoch=$epoch (best epoch=$best_epoch val_acc=$(round(best_val_acc; digits=3)))"
                    break
                end
            end
        elseif config.verbose
            @info "epoch=$epoch train=$(round(epoch_loss/n_train; digits=4))/$(round(epoch_correct/n_train; digits=3))"
        end
    end

    if config.restore_best && best_state !== nothing
        Flux.loadmodel!(model, best_state)
        config.verbose && @info "restored best checkpoint" best_epoch best_val_acc=round(best_val_acc; digits=3)
    end
    return history
end
