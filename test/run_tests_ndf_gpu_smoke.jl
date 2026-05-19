#
# GPU smoke test for NDF. Verifies that the model + Â + per-batch tensors
# can live on GPU, that forward/backward complete without scalar-indexing
# errors, and that the GPU forward matches the CPU forward to ≤1e-4.
#
# Skip with NDF_GPU_SKIP=1 (used in CI without a GPU).
#
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Random, SparseArrays, LinearAlgebra
using Flux
using CUDA, cuDNN
using Test
using Adjacently
using Adjacently.Fingerprints: NDF, prepare_adjacency, default_node_features

if !CUDA.functional()
    @info "CUDA not functional — skipping GPU smoke test"
    exit(0)
end

CUDA.allowscalar(false)
@info "GPU smoke test" device=CUDA.name(CUDA.device()) free_GB=round(Int, CUDA.available_memory() / 2^30)

# Build a small synthetic graph + features deterministically.
Random.seed!(0)
n = 200
density = 0.02
nnz = Int(round(n * n * density))
rows = rand(1:n, nnz); cols = rand(1:n, nnz); vals = ones(Float32, nnz)
A = sparse(rows, cols, vals, n, n)
A = max.(A, A')                 # symmetrize
foreach(i -> A[i,i] = 0, 1:n)   # zero diag (prepare_adjacency adds self-loops)
dropzeros!(A)

Â_cpu = prepare_adjacency(A)
X_cpu = default_node_features(A)
n_vertices = n
d_in = 1 + size(X_cpu, 2)

# Build a tiny NDF identical on CPU and GPU.
Random.seed!(42)
ndf_cpu = NDF(d_in, 2; hidden=2, K=5, α=0.15f0, dropout=0.0f0,
              readout=:flatten, n_vertices=n_vertices)
Flux.testmode!(ndf_cpu)

# Forward on CPU with a batch of 4 docs.
B = 4
Φ_cpu = zeros(Float32, n, d_in, B)
mask_cpu = zeros(Bool, n, B)
for b in 1:B
    seeds = rand(1:n, 5)
    for u in seeds
        Φ_cpu[u, 1, b] = 1.0f0 / 5
        mask_cpu[u, b] = true
    end
    @views Φ_cpu[:, 2:end, b] .= X_cpu
end
out_cpu = ndf_cpu(Φ_cpu, Â_cpu; seed_mask=mask_cpu)
@info "CPU forward OK" size_out=size(out_cpu)

# Move everything to GPU and re-run.
ndf_gpu = Flux.gpu(ndf_cpu)
Â_gpu = Flux.gpu(Â_cpu)
X_gpu = Flux.gpu(X_cpu)
Φ_gpu = Flux.gpu(Φ_cpu)
mask_gpu = Flux.gpu(mask_cpu)
@info "Tensors on GPU" Â_type=string(typeof(Â_gpu)) Φ_type=string(typeof(Φ_gpu))

out_gpu = Array(ndf_gpu(Φ_gpu, Â_gpu; seed_mask=mask_gpu))
@info "GPU forward OK" size_out=size(out_gpu)

max_diff = maximum(abs.(out_cpu .- out_gpu))
@info "CPU vs GPU forward agreement" max_abs_diff=max_diff

@testset "GPU smoke" begin
    @test size(out_gpu) == size(out_cpu)
    @test max_diff < 1f-4
end

# Backward smoke: grad of a scalar loss w.r.t. GPU model params.
Flux.trainmode!(ndf_gpu)
y_gpu = Flux.gpu(Flux.onehotbatch([1, 2, 1, 2], 1:2))
loss, grads = Flux.withgradient(ndf_gpu) do m
    Flux.logitcrossentropy(m(Φ_gpu, Â_gpu; seed_mask=mask_gpu), y_gpu)
end
@info "GPU backward OK" loss=Float32(loss)

@testset "GPU backward" begin
    @test isfinite(Float32(loss))
end

@info "GPU smoke test complete"
