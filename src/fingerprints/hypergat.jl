#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2026 Jimmy Dubuisson <jimmy@dubuisson.ch>
#
# GNU GPL v3 or later.
#

# =============================================================================
# HyperGAT (Ding et al., EMNLP 2020) — per-document hypergraph with two-level
# attention (node→hyperedge, hyperedge→node). The clique-expansion shortcut
# (hypergraph.jl) was null on Ohsumed because it discards higher-order structure;
# this is the faithful mechanism (published Ohsumed 69.82, one-hot, no pretrain).
#
# Hyperedges per document:
#   - sequential: fixed-size chunks of the token stream (sentence proxy),
#   - semantic:   one per LDA topic present in the doc (its top-word members).
#
# Everything is expressed as (sparse incidence × dense) products so it runs on
# the GPU through the same materialized-tangent rrule used by NDF. Rectangular
# incidence means the backward needs an EXPLICIT transpose operator (lazy
# Transpose{CuSparseCSR}*dense has no GPU method) — hence `_hspmm(A, At, X)`
# carries both directions as real CSR matrices.
# =============================================================================

# A*X forward; backward uses the explicitly-provided transpose At (both real CSR).
_hspmm(A::AbstractMatrix, At::AbstractMatrix, X::AbstractMatrix) = A * X
function ChainRulesCore.rrule(::typeof(_hspmm), A::AbstractMatrix,
                              At::AbstractMatrix, X::AbstractMatrix)
    Y = A * X
    function _hspmm_pb(ΔY)
        ΔYd = _spmm_tangent_like(ΔY, Y)          # materialize on Y's backend
        return (NoTangent(), NoTangent(), NoTangent(), At * ΔYd)
    end
    return Y, _hspmm_pb
end

# ── Per-document hypergraph ──────────────────────────────────────────────────
"""
    HyperDoc

One document as a hypergraph. `node_ids` = distinct global vocab ids; incidence
is the pair list `(pair_node, pair_edge)` in LOCAL indices (node 1..L, edge
1..M); `n_edges` = M; `seed_w` per-node lexical scalar; `label`.
"""
struct HyperDoc
    node_ids::Vector{Int32}
    pair_node::Vector{Int32}
    pair_edge::Vector{Int32}
    n_edges::Int
    seed_w::Vector{Float32}
    label::Int
end

"""
    build_hyper_doc(tokens, vocab, topic_sets; seq_chunk=20, weight=:freq) -> HyperDoc | nothing

Sequential hyperedges = consecutive `seq_chunk`-token chunks; semantic hyperedges
= one per topic whose top-words intersect the doc (≥2 members). Returns `nothing`
for docs with no in-vocab tokens.
"""
function build_hyper_doc(tokens, vocab::Dict{String,Int},
                         topic_sets::Vector{Set{Int32}};
                         seq_chunk::Int=20, weight::Symbol=:freq)
    seq = Int32[]
    for t in tokens
        id = get(vocab, String(t), 0); id == 0 && continue
        push!(seq, Int32(id))
    end
    isempty(seq) && return nothing

    local_of = Dict{Int32,Int32}(); node_ids = Int32[]; counts = Float32[]
    for g in seq
        li = get(local_of, g, Int32(0))
        if li == 0
            push!(node_ids, g); push!(counts, 1f0); local_of[g] = Int32(length(node_ids))
        else
            counts[li] += 1f0
        end
    end
    L = length(node_ids)

    pair_node = Int32[]; pair_edge = Int32[]; eidx = 0
    # sequential hyperedges (distinct nodes within each chunk)
    n = length(seq)
    p = 1
    while p <= n
        q = min(n, p + seq_chunk - 1)
        members = Set{Int32}()
        for r in p:q; push!(members, local_of[seq[r]]); end
        if length(members) >= 1
            eidx += 1
            for m in members; push!(pair_node, m); push!(pair_edge, Int32(eidx)); end
        end
        p = q + 1
    end
    # semantic hyperedges (doc nodes that are top-words of a topic)
    for S in topic_sets
        members = Int32[]
        @inbounds for i in 1:L; node_ids[i] in S && push!(members, Int32(i)); end
        if length(members) >= 2
            eidx += 1
            for m in members; push!(pair_node, m); push!(pair_edge, Int32(eidx)); end
        end
    end
    eidx == 0 && return nothing

    total = sum(counts)
    seed_w = weight === :binary ? ones(Float32, L) : (counts ./ total)
    return HyperDoc(node_ids, pair_node, pair_edge, eidx, seed_w, 0)
end

function docs_to_hyper_docs(docs::Vector{TextGcnDoc}, vocab::Dict{String,Int},
                            label2id::Dict{String,Int}, topics::Vector{Vector{Int32}};
                            seq_chunk::Int=20, weight::Symbol=:freq)
    topic_sets = [Set(t) for t in topics if length(t) >= 2]
    train = HyperDoc[]; test = HyperDoc[]
    for d in docs
        hd = build_hyper_doc(d.tokens, vocab, topic_sets; seq_chunk=seq_chunk, weight=weight)
        hd === nothing && continue
        hd = HyperDoc(hd.node_ids, hd.pair_node, hd.pair_edge, hd.n_edges, hd.seed_w, label2id[d.label])
        (d.split === :train ? train : test) |> v -> push!(v, hd)
    end
    return train, test, length(label2id)
end

# ── Batch (block-diagonal incidence) ─────────────────────────────────────────
struct HyperBatch{NI,S,SW,SM,SZ}
    node_ids::NI     # (N,)
    V2P::S; V2Pt::S  # (N×P),(P×N) node–pair incidence + transpose
    E2P::S; E2Pt::S  # (M×P),(P×M) edge–pair incidence + transpose
    seed_w::SW       # (N,)
    Ind::SM; Indt::SM # (N×B),(B×N) node→doc segment
    sizes::SZ        # (B,)
    labels::Vector{Int}
end

function build_hyper_batch(hds::AbstractVector{HyperDoc})
    B = length(hds)
    Ls = Int[length(h.node_ids) for h in hds]
    Ms = Int[h.n_edges for h in hds]
    noff = cumsum(vcat(0, Ls[1:end-1])); eoff = cumsum(vcat(0, Ms[1:end-1]))
    N = sum(Ls); M = sum(Ms)
    node_ids = Vector{Int32}(undef, N); seed_w = Vector{Float32}(undef, N)
    vI = Int[]; vJ = Int[]; eI = Int[]; eJ = Int[]         # incidence COO (global pair idx)
    indR = Int[]; indC = Int[]; labels = Vector{Int}(undef, B)
    poff = 0
    for (b, h) in enumerate(hds)
        no = noff[b]; eo = eoff[b]; L = Ls[b]
        @inbounds for i in 1:L
            node_ids[no+i] = h.node_ids[i]; seed_w[no+i] = h.seed_w[i]
            push!(indR, no+i); push!(indC, b)
        end
        @inbounds for k in eachindex(h.pair_node)
            gp = poff + k
            push!(vI, no + Int(h.pair_node[k])); push!(vJ, gp)
            push!(eI, eo + Int(h.pair_edge[k])); push!(eJ, gp)
        end
        poff += length(h.pair_node)
    end
    P = poff
    V2P = sparse(vI, vJ, ones(Float32, P), N, P)
    E2P = sparse(eI, eJ, ones(Float32, P), M, P)
    Ind = sparse(indR, indC, ones(Float32, N), N, B)
    HyperBatch(node_ids, V2P, sparse(V2P'), E2P, sparse(E2P'),
               seed_w, Ind, sparse(Ind'), Float32.(Ls), labels_from(hds))
end
labels_from(hds) = Int[h.label for h in hds]

# ── Model ────────────────────────────────────────────────────────────────────
struct HyperGATLayer{M,A}
    W1::M; W2::M; a1::A; a2::A
end
Flux.@layer HyperGATLayer
function HyperGATLayer(d_in::Int, d_out::Int)
    g(a,b) = Float32.((rand(Float32,a,b) .- 0.5f0) .* (2f0*sqrt(6f0/(a+b))))
    HyperGATLayer(g(d_in,d_out), g(d_out,d_out), g(d_out,1), g(d_out,1))
end

_lrelu(x) = max.(x, 0.2f0 .* x)

# X: (N×d_in). Returns (N×d_out). b = (V2P,V2Pt,E2P,E2Pt).
function (L::HyperGATLayer)(X, V2P, V2Pt, E2P, E2Pt)
    WH  = X * L.W1                              # N×d1
    WHp = _hspmm(V2Pt, V2P, WH)                 # P×d1  (gather node→pair)
    s1  = vec(_lrelu(WHp) * L.a1)               # P
    e1  = exp.(s1 .- maximum(s1))               # P
    esum = _hspmm(E2P, E2Pt, reshape(e1, :, 1)) # M×1  (sum per edge)
    den1 = vec(_hspmm(E2Pt, E2P, esum))         # P    (back to pairs)
    α   = e1 ./ (den1 .+ 1f-9)
    F   = relu.(_hspmm(E2P, E2Pt, α .* WHp))    # M×d1 (edge repr)
    W2F  = F * L.W2                             # M×d_out
    W2Fp = _hspmm(E2Pt, E2P, W2F)               # P×d_out (gather edge→pair)
    s2  = vec(_lrelu(W2Fp) * L.a2)              # P
    e2  = exp.(s2 .- maximum(s2))
    nsum = _hspmm(V2P, V2Pt, reshape(e2, :, 1)) # N×1 (sum per node)
    den2 = vec(_hspmm(V2Pt, V2P, nsum))         # P
    β   = e2 ./ (den2 .+ 1f-9)
    relu.(_hspmm(V2P, V2Pt, β .* W2Fp))         # N×d_out (node repr)
end

struct HyperGAT{EM,L1,L2,DR,CL}
    emb::EM; l1::L1; l2::L2; drop::DR; classifier::CL
end
Flux.@layer HyperGAT trainable=(emb, l1, l2, classifier)

function HyperGAT(V::Int, n_classes::Int; d_emb::Int=200, h1::Int=200, h2::Int=100,
                  dropout::Float32=0.5f0, head_hidden::Int=256, emb_init=nothing)
    g(a,b) = Float32.((rand(Float32,a,b) .- 0.5f0) .* (2f0*sqrt(6f0/(a+b))))
    emb = emb_init === nothing ? g(d_emb, V) : Float32.(emb_init)
    cls = Chain(Dense(h2 => head_hidden, relu), Flux.Dropout(dropout), Dense(head_hidden => n_classes))
    HyperGAT(emb, HyperGATLayer(d_emb, h1), HyperGATLayer(h1, h2), Flux.Dropout(dropout), cls)
end

function (m::HyperGAT)(b::HyperBatch)
    X = m.emb[:, b.node_ids]'                       # N×d_emb
    H = m.drop(m.l1(X, b.V2P, b.V2Pt, b.E2P, b.E2Pt))
    H = m.drop(m.l2(H, b.V2P, b.V2Pt, b.E2P, b.E2Pt))
    pooled = (b.Indt * H) ./ b.sizes                # B×h2 (mean over doc nodes)
    m.classifier(pooled')                           # n_classes×B
end

# ── Training (mirrors train_perdoc!) ─────────────────────────────────────────
function train_hypergat!(model::HyperGAT, train_docs::AbstractVector{HyperDoc};
                         val_docs=nothing, test_docs=nothing, n_classes::Int,
                         config::TrainConfig=TrainConfig(),
                         movebatch::Function=identity, movearr::Function=identity,
                         freeze_emb::Bool=false)
    opt = config.weight_decay > 0f0 ?
        Flux.AdamW(config.lr, (0.9f0,0.999f0), config.weight_decay) : Flux.Adam(config.lr)
    st = Flux.setup(opt, model)
    freeze_emb && Flux.Optimisers.freeze!(st.emb)
    rng = MersenneTwister(config.seed)
    n = length(train_docs)
    hist = (train_loss=Float32[], train_acc=Float32[], val_acc=Float32[], test_acc=Float32[])
    best = config.initial_best_val_acc; best_state = nothing; bepoch = 0; pat = 0
    for epoch in 1:config.epochs
        Flux.trainmode!(model)
        perm = Random.shuffle(rng, collect(1:n)); eloss = 0f0; ecorr = 0
        for s in 1:config.batch_size:n
            e = min(s+config.batch_size-1, n)
            cb = build_hyper_batch(train_docs[perm[s:e]]); b = movebatch(cb)
            y = movearr(Flux.onehotbatch(cb.labels, 1:n_classes))
            loss, gs = Flux.withgradient(m -> Flux.logitcrossentropy(m(b), y), model)
            Flux.update!(st, model, gs[1])
            eloss += Float32(loss)*length(cb.labels)
            ŷ = Array(model(b)); ecorr += count(vec(map(I->I[1], argmax(ŷ;dims=1))) .== cb.labels)
        end
        push!(hist.train_loss, eloss/n); push!(hist.train_acc, ecorr/n)
        if val_docs !== nothing && epoch % config.val_every == 0
            va = _eval_hypergat(model, val_docs, config.batch_size, n_classes, movebatch)
            push!(hist.val_acc, va)
            config.verbose && @info "epoch=$epoch train=$(round(eloss/n;digits=4))/$(round(ecorr/n;digits=3)) val=$(round(va;digits=4))"
            if va > best + 1f-5
                best = va; bepoch = epoch; config.restore_best && (best_state = Flux.state(model)); pat = 0
                if test_docs !== nothing
                    ta = _eval_hypergat(model, test_docs, config.batch_size, n_classes, movebatch)
                    push!(hist.test_acc, ta); config.verbose && @info "  test epoch=$epoch test_acc=$(round(ta;digits=4))"
                end
            else
                test_docs !== nothing && push!(hist.test_acc, NaN32); pat += 1
                config.early_stop_patience > 0 && pat >= config.early_stop_patience && (@info "early stop epoch=$epoch"; break)
            end
        end
    end
    if config.restore_best && best_state !== nothing
        Flux.loadmodel!(model, best_state); @info "restored best" bepoch best_val=round(best;digits=3)
    end
    return hist
end

function _eval_hypergat(model, docs, bs, n_classes, movebatch)
    Flux.testmode!(model); correct = 0; n = length(docs)
    for s in 1:bs:n
        e = min(s+bs-1, n); cb = build_hyper_batch(docs[s:e]); b = movebatch(cb)
        ŷ = Array(model(b)); correct += count(vec(map(I->I[1], argmax(ŷ;dims=1))) .== cb.labels)
    end
    Flux.trainmode!(model); correct / n
end
