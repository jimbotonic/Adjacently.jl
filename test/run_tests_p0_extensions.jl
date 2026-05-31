#
# Smoke test for the P0 extensions from GPT_SUGGESTIONS.md:
#   - weighted PPMI graph (build_domain_graph; weighted=true)
#   - TF-IDF seed weighting (textgcn_to_documents; weight=:tfidf)
#   - OPC central-node projection (NDF readout=:opc_flatten + top_central_nodes)
#
# Runs on a tiny synthetic corpus so it stays CPU-fast (<5 s).
#
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Random, SparseArrays, LinearAlgebra
using Flux
using Test
using Adjacently
using Adjacently.Fingerprints:
    NDF, Document, TrainConfig, train!,
    prepare_adjacency, default_node_features, top_central_nodes,
    build_domain_graph,
    TextGcnDoc, textgcn_to_documents

Random.seed!(0)

# --- 1. weighted PPMI graph preserves magnitudes ----------------------------
n = 50
K = sprand(Float32, n, n, 0.05)  # random PPMI-like values in [0, 1]
K = max.(K, K')                  # symmetrize
foreach(i -> K[i, i] = 0f0, 1:n)
dropzeros!(K)

A_bin = build_domain_graph(K, 0.02)
A_wt  = build_domain_graph(K, 0.02; weighted=true)

@testset "build_domain_graph weighted flag" begin
    @test size(A_bin) == size(A_wt)
    @test length(A_bin.nzval) == length(A_wt.nzval)
    @test all(A_bin.nzval .== 1.0f0)
    @test any(A_wt.nzval .!= 1.0f0)        # at least one non-unit weight
    @test all(A_wt.nzval .> 0.0f0)         # all positive (PPMI ≥ 0)
    @test maximum(A_wt.nzval) <= maximum(K.nzval) + 1f-6   # bounded by source
end

# --- 2. TF-IDF seed weighting -----------------------------------------------
# Tiny corpus: 4 train docs, 2 classes, 1 test doc.
docs = [
    TextGcnDoc(["the", "cat", "sat"],         "a", :train),
    TextGcnDoc(["the", "dog", "ran"],         "b", :train),
    TextGcnDoc(["the", "cat", "ran"],         "a", :train),
    TextGcnDoc(["the", "dog", "sat"],         "b", :train),
    TextGcnDoc(["the", "cat", "ran"],         "a", :test),
]
vocab = Dict("the" => 1, "cat" => 2, "dog" => 3, "sat" => 4, "ran" => 5)
label2id = Dict("a" => 1, "b" => 2)

tr_freq,  _, _ = textgcn_to_documents(docs, vocab, label2id; weight=:freq)
tr_tfidf, _, _ = textgcn_to_documents(docs, vocab, label2id; weight=:tfidf)

@testset "TF-IDF seed weighting" begin
    @test length(tr_freq) == length(tr_tfidf) == 4
    # In doc 1 ["the","cat","sat"], "the" should get LESS weight under tfidf
    # than freq because it appears in all 4 training docs (df=4) while
    # "cat" appears in only 2 (df=2).
    freq_doc1  = Dict(tr_freq[1].seed_indices  .=> tr_freq[1].seed_weights)
    tfidf_doc1 = Dict(tr_tfidf[1].seed_indices .=> tr_tfidf[1].seed_weights)
    the_id, cat_id = vocab["the"], vocab["cat"]
    @test freq_doc1[the_id] ≈ freq_doc1[cat_id]            # equal counts → equal freq
    @test tfidf_doc1[the_id] < tfidf_doc1[cat_id]          # tfidf downweights common
    # Both representations sum to 1 (probability distribution).
    @test sum(tr_freq[1].seed_weights) ≈ 1f0 atol=1f-5
    @test sum(tr_tfidf[1].seed_weights) ≈ 1f0 atol=1f-5
end

# --- 3. OPC central-node projection -----------------------------------------
# Build a small symmetric graph; top_central_nodes should rank vertices.
n2 = 60
A2 = sprand(Float32, n2, n2, 0.1)
A2 = max.(A2, A2')
foreach(i -> A2[i, i] = 0f0, 1:n2)
dropzeros!(A2)
Â = prepare_adjacency(A2)

central = top_central_nodes(Â, 8)
@testset "top_central_nodes" begin
    @test length(central) == 8
    @test length(unique(central)) == 8       # no duplicates
    @test all(1 .<= central .<= n2)
end

# Build an NDF with :opc_flatten readout and run forward + backward on a
# tiny batch. Verifies the new readout path doesn't break the existing
# encoder + propagate + classifier chain.
d_in = 1 + 4                                 # seed weight + 4 degree features
hidden = 2
B = 3
ndf = NDF(d_in, 2; hidden=hidden, K=5, α=0.15f0, dropout=0.0f0,
          readout=:opc_flatten, n_central=length(central))
Flux.trainmode!(ndf)

X = default_node_features(A2)                # (n2 × 4)
Φ = zeros(Float32, n2, d_in, B)
mask = zeros(Bool, n2, B)
for b in 1:B
    seeds = rand(1:n2, 4)
    for u in seeds
        Φ[u, 1, b] = 1.0f0 / 4
        mask[u, b] = true
    end
    @views Φ[:, 2:end, b] .= X
end
y = Flux.onehotbatch([1, 2, 1], 1:2)

out = ndf(Φ, Â; seed_mask=mask, central_nodes=central)
@testset "NDF :opc_flatten forward" begin
    @test size(out) == (2, B)                # n_classes × B
    @test all(isfinite, out)
end

loss, grads = Flux.withgradient(ndf) do m
    Flux.logitcrossentropy(m(Φ, Â; seed_mask=mask, central_nodes=central), y)
end
@testset "NDF :opc_flatten backward" begin
    @test isfinite(Float32(loss))
    # Classifier weights should be (2 × n_central*hidden) = (2 × 16).
    @test size(ndf.classifier[2].weight) == (2, length(central) * hidden)
end

@info "P0 extensions smoke test complete"
