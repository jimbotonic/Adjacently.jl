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
end

function TrainConfig(;
        proxy_epochs::Int=100,
        proxy_lr::Float64=0.001,
        reinforce_epochs::Int=50,
        reinforce_lr::Float64=0.0001,
        sigma::Float64=0.1,
        baseline_ema::Float64=0.9)
    return TrainConfig(proxy_epochs, proxy_lr, reinforce_epochs, reinforce_lr, sigma, baseline_ema)
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
        raw_scores = gnn_forward(model)

        # Normalize to zero mean, unit variance
        mu = sum(raw_scores) / n
        centered = raw_scores .- mu
        sigma_sq = sum(centered .^ 2) / n
        if sigma_sq < 1e-12
            raw_scores .+= randn(n) .* 0.01
            mu = sum(raw_scores) / n
            centered = raw_scores .- mu
            sigma_sq = sum(centered .^ 2) / n
        end
        sigma_val = sqrt(sigma_sq)
        scores = centered ./ sigma_val

        # Proxy loss and gradient
        loss = 0.0
        dL_ds = zeros(n)
        for v in vertices(g)
            vi = Int(v)
            sv = scores[vi]
            for u in outneighbors(g, v)
                ui = Int(u)
                diff = sv - scores[ui]
                loss += diff^2
                dL_ds[vi] += 2.0 * diff
                dL_ds[ui] -= 2.0 * diff
            end
        end

        mean_dL = sum(dL_ds) / n
        score_grad_dot = sum(scores .* dL_ds) / n
        dL_dZ2 = (dL_ds .- mean_dL) ./ sigma_val .- scores .* score_grad_dot ./ sigma_val

        push!(losses, loss)
        if epoch == 1 || epoch % 10 == 0
            @info "Proxy epoch $epoch/$(config.proxy_epochs): loss = $(round(loss, digits=4)), score_std = $(round(sigma_val, digits=6))"
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

