#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2026 Jimmy Dubuisson <jimmy@dubuisson.ch>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#

# =============================================================================
# HyperGAT-style semantic (LDA-topic) hyperedges for the per-document
# graph-of-words.
#
# The per-doc pairwise/window graph (perdoc.jl) captures *sequential* within-doc
# structure; the Phase-0 controls showed that alone ties/loses to bag-of-words on
# topical data. HyperGAT (Ding et al., EMNLP 2020) beats BoW without pretraining
# (Ohsumed 69.82, one-hot) by adding *semantic hyperedges*: an LDA topic is a
# hyperedge connecting the document's words that belong to that topic — genuine
# higher-order structure that pairwise co-occurrence cannot express.
#
# We fold that signal into the existing pipeline by CLIQUE-EXPANDING each semantic
# hyperedge into word-word edges, so `PerDocGNN` / `train_perdoc!` / the GPU-safe
# `_spmm` batching all consume it unchanged. LDA is fit on the TRAIN split only
# (inductive); test documents reuse the fitted topics.
# =============================================================================

"""
    lda_topic_words(train_docs, vocab; K, iters, top_n, α, β, seed)
        -> Vector{Vector{Int32}}

Collapsed-Gibbs LDA over the TRAIN split (token sequences of in-vocab word ids).
Returns, for each of `K` topics, its `top_n` highest-count word ids — the member
set of that topic's semantic hyperedge. Train-only ⇒ inductive.
"""
function lda_topic_words(train_docs::Vector{TextGcnDoc}, vocab::Dict{String,Int};
                         K::Int=50, iters::Int=150, top_n::Int=15,
                         α::Float32=0.1f0, β::Float32=0.01f0, seed::Int=0)
    V = length(vocab)
    seqs = Vector{Vector{Int32}}()
    for d in train_docs
        d.split === :train || continue
        s = Int32[]
        for t in d.tokens
            id = get(vocab, String(t), 0); id == 0 && continue
            push!(s, Int32(id))
        end
        isempty(s) || push!(seqs, s)
    end
    D = length(seqs)
    D == 0 && return [Int32[] for _ in 1:K]

    rng = MersenneTwister(seed)
    n_dk = zeros(Int32, D, K)
    n_kw = zeros(Int32, K, V)
    n_k  = zeros(Int32, K)
    z = [Vector{Int32}(undef, length(s)) for s in seqs]
    @inbounds for d in 1:D, i in eachindex(seqs[d])
        k = rand(rng, 1:K); z[d][i] = Int32(k); w = seqs[d][i]
        n_dk[d, k] += 1; n_kw[k, w] += 1; n_k[k] += 1
    end

    Vβ = Float32(V) * β
    p = zeros(Float32, K)
    for _ in 1:iters
        @inbounds for d in 1:D
            sd = seqs[d]; zd = z[d]
            for i in eachindex(sd)
                w = sd[i]; k0 = zd[i]
                n_dk[d, k0] -= 1; n_kw[k0, w] -= 1; n_k[k0] -= 1
                s = 0f0
                for k in 1:K
                    p[k] = (n_dk[d, k] + α) * (n_kw[k, w] + β) / (n_k[k] + Vβ)
                    s += p[k]
                end
                r = rand(rng) * s; acc = 0f0; knew = K
                for k in 1:K
                    acc += p[k]
                    if r <= acc; knew = k; break; end
                end
                zd[i] = Int32(knew)
                n_dk[d, knew] += 1; n_kw[knew, w] += 1; n_k[knew] += 1
            end
        end
    end

    topics = Vector{Vector{Int32}}(undef, K)
    for k in 1:K
        row = collect(@view n_kw[k, :])
        m = min(top_n, V)
        idx = partialsortperm(row, 1:m; rev=true)
        topics[k] = Int32.(idx)
    end
    return topics
end

# Clique-expand the doc's semantic hyperedges into its edge set. `topic_sets` are
# the topics' top-word id sets; for each topic, the doc's member nodes form a
# clique with weight `sem_w` (accumulated onto any existing window edge).
function _add_semantic_edges!(ew::Dict{Tuple{Int32,Int32},Float32},
                              node_ids::Vector{Int32},
                              topic_sets::Vector{Set{Int32}}, sem_w::Float32)
    loc = Dict{Int32,Int32}()
    @inbounds for i in eachindex(node_ids); loc[node_ids[i]] = Int32(i); end
    for S in topic_sets
        members = Int32[]
        @inbounds for i in eachindex(node_ids)
            node_ids[i] in S && push!(members, Int32(i))
        end
        length(members) < 2 && continue
        @inbounds for ii in 1:length(members)-1, jj in ii+1:length(members)
            a = min(members[ii], members[jj]); b = max(members[ii], members[jj])
            ew[(a, b)] = get(ew, (a, b), 0f0) + sem_w
        end
    end
end

"""
    augment_with_semantics(lg, topic_sets, sem_w) -> LocalGraph

Add clique-expanded semantic (LDA-topic) edges to a window `LocalGraph`.
"""
function augment_with_semantics(lg::LocalGraph, topic_sets::Vector{Set{Int32}},
                                sem_w::Float32)
    ew = Dict{Tuple{Int32,Int32},Float32}()
    @inbounds for k in eachindex(lg.src)
        a = min(lg.src[k], lg.dst[k]); b = max(lg.src[k], lg.dst[k])
        ew[(a, b)] = get(ew, (a, b), 0f0) + lg.ew[k]
    end
    _add_semantic_edges!(ew, lg.node_ids, topic_sets, sem_w)
    src = Int32[]; dst = Int32[]; w = Float32[]
    for ((a, b), c) in ew
        push!(src, a); push!(dst, b); push!(w, c)
    end
    return LocalGraph(lg.node_ids, src, dst, w, lg.seed_w, lg.label)
end

"""
    docs_to_hyper_graphs(docs, vocab, label2id, topics; window, weight, sem_w)
        -> (train, test, n_classes)

Like `docs_to_local_graphs`, but each per-doc window graph is augmented with
clique-expanded semantic hyperedges from the fitted LDA `topics` (a
`Vector{Vector{Int32}}` of top-word ids per topic).
"""
function docs_to_hyper_graphs(docs::Vector{TextGcnDoc},
                              vocab::Dict{String,Int},
                              label2id::Dict{String,Int},
                              topics::Vector{Vector{Int32}};
                              window::Int=3, weight::Symbol=:freq,
                              sem_w::Float32=1f0)
    topic_sets = [Set(t) for t in topics if length(t) >= 2]
    train = LocalGraph[]; test = LocalGraph[]
    for d in docs
        lg = build_local_window_graph(d.tokens, vocab; window=window, weight=weight)
        lg === nothing && continue
        lg = LocalGraph(lg.node_ids, lg.src, lg.dst, lg.ew, lg.seed_w, label2id[d.label])
        lg = augment_with_semantics(lg, topic_sets, sem_w)
        (d.split === :train ? train : test) |> v -> push!(v, lg)
    end
    return train, test, length(label2id)
end
