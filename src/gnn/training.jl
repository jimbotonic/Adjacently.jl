using LinearAlgebra: norm

"""
    TrainConfig

Configuration for GNN training phases.
"""
struct TrainConfig
    proxy_epochs::Int
    proxy_lr::Float64
    reinforce_epochs::Int
    reinforce_lr::Float64
    sigma::Float64          # Exploration noise std for REINFORCE
    baseline_ema::Float64   # EMA decay for baseline (0.9 typical)
    grad_clip_norm::Float64 # Max norm for gradient clipping
end

function TrainConfig(;
        proxy_epochs::Int=100,
        proxy_lr::Float64=0.001,
        reinforce_epochs::Int=50,
        reinforce_lr::Float64=0.0001,
        sigma::Float64=0.1,
        baseline_ema::Float64=0.9,
        grad_clip_norm::Float64=1.0)
    return TrainConfig(proxy_epochs, proxy_lr, reinforce_epochs, reinforce_lr, sigma, baseline_ema, grad_clip_norm)
end

"""
    train_gnn_proxy!(model::GNNModel, g::AbstractGraph{T}, config::TrainConfig) where {T<:Unsigned}

Phase 1: Train GNN with spectral ordering proxy loss.
Minimizes bandwidth L = Σ_{(u,v) ∈ E} (s[u] - s[v])² subject to unit variance.
"""
function train_gnn_proxy!(model::GNNModel, g::AbstractGraph{T}, config::TrainConfig) where {T<:Unsigned}
    losses = Float64[]
    n = model.n

    for epoch in 1:config.proxy_epochs
        raw_scores = gnn_forward(model; training=false)

        # Raw score statistics
        mu = sum(raw_scores) / n
        centered = raw_scores .- mu
        sigma_sq = sum(centered .^ 2) / n
        sigma_val = sqrt(max(sigma_sq, 1e-12))

        # Symmetric Laplacian action on raw scores and bandwidth
        bandwidth = 0.0
        Ls = zeros(n)
        for v in vertices(g)
            vi = Int(v)
            for u in outneighbors(g, v)
                ui = Int(u)
                diff = raw_scores[vi] - raw_scores[ui]
                bandwidth += diff^2
                Ls[vi] += diff
                Ls[ui] -= diff
            end
        end

        # Rayleigh quotient: R = bandwidth / (n · σ²)
        # Minimizing R converges to the Fiedler vector.
        # Gradient dR ∝ Ls - R·centered: bandwidth pulls neighbors together,
        # R·centered pushes scores apart. At the Fiedler eigenvector these
        # balance and the gradient is zero. Variance is preserved (first order)
        # because dot(centered, Ls - R·centered) = bandwidth - R·n·σ² = 0.
        rayleigh = bandwidth / (n * max(sigma_sq, 1e-12))
        dL_dZ2 = Ls .- rayleigh .* centered

        # Gradient clipping
        grad_norm = norm(dL_dZ2)
        if grad_norm > config.grad_clip_norm
            dL_dZ2 .*= (config.grad_clip_norm / grad_norm)
        end

        push!(losses, rayleigh)
        if epoch <= 5 || epoch % 10 == 0
            @info "Proxy epoch $epoch/$(config.proxy_epochs): loss = $(round(rayleigh, digits=4)), score_std = $(round(sigma_val, digits=6))"
        end

        gnn_backward!(model, dL_dZ2; lr=config.proxy_lr)
    end

    return losses
end

"""
    train_gnn_reinforce!(model::GNNModel, g::AbstractGraph{T}, compress_fn,
                         config::TrainConfig) where {T<:Unsigned}

Phase 2: REINFORCE training using actual compression bpe as reward signal.
`compress_fn(ordering::Dict{T,T})::Float64` is a closure returning bits-per-edge.
"""
function train_gnn_reinforce!(model::GNNModel, g::AbstractGraph{T}, compress_fn,
                              config::TrainConfig) where {T<:Unsigned}
    results = Tuple{Int,Float64}[]
    baseline_bpe = nothing

    for epoch in 1:config.reinforce_epochs
        scores = gnn_forward(model)
        noise = randn(model.n) .* config.sigma
        noisy_scores = scores .+ noise

        perm = sortperm(noisy_scores)
        ordering = Dict(T(old) => T(new) for (new, old) in enumerate(perm))
        bpe = compress_fn(ordering)

        if baseline_bpe === nothing
            baseline_bpe = bpe
        else
            baseline_bpe = config.baseline_ema * baseline_bpe + (1.0 - config.baseline_ema) * bpe
        end

        reward = -(bpe - baseline_bpe)
        dL_dZ2 = -reward .* noise ./ (config.sigma^2)

        push!(results, (epoch, bpe))
        @info "REINFORCE epoch $epoch/$(config.reinforce_epochs): bpe = $(round(bpe, digits=4)), baseline = $(round(baseline_bpe, digits=4))"
        gnn_backward!(model, dL_dZ2; lr=config.reinforce_lr)
    end

    return results
end

