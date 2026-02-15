#
# GNN-augmented actor-critic: linear heads over per-vertex GNN features.
# The features are provided by the caller via a function feature_fn(v)::Vector{Float64}.
#

using LinearAlgebra: dot

mutable struct GnnACPolicy
    Wpi::Matrix{Float64}   # (feat_dim × num_actions)
    bpi::Vector{Float64}   # (num_actions)
    Wv::Vector{Float64}    # (feat_dim)
    bv::Float64
    actor_lr::Float64
    critic_lr::Float64
    gamma::Float64
    temperature::Float64
    feat_dim::Int
    num_actions::Int
end

function GnnACPolicy(feat_dim::Int; num_actions::Int=NUM_ACTIONS,
                     actor_lr::Float64=0.05, critic_lr::Float64=0.1,
                     gamma::Float64=0.0, temperature::Float64=1.0)
    Wpi = zeros(Float64, feat_dim, num_actions)
    bpi = zeros(Float64, num_actions)
    Wv = zeros(Float64, feat_dim)
    bv = 0.0
    return GnnACPolicy(Wpi, bpi, Wv, bv, actor_lr, critic_lr, gamma, temperature, feat_dim, num_actions)
end

"""
    save_gnn_ac_policy(policy, filepath)

Save GnnACPolicy weights (Wpi, bpi, Wv, bv) in binary format.
"""
function save_gnn_ac_policy(policy::GnnACPolicy, filepath::AbstractString)
    open(filepath, "w") do f
        write(f, Int32(policy.feat_dim))
        write(f, Int32(policy.num_actions))
        write(f, policy.Wpi)
        write(f, policy.bpi)
        write(f, policy.Wv)
        write(f, Float64(policy.bv))
    end
end

"""
    load_gnn_ac_policy(filepath; actor_lr=0.05, critic_lr=0.1, gamma=0.0, temperature=1.0)

Load GnnACPolicy weights from file. Learning rate and temperature are set from kwargs.
"""
function load_gnn_ac_policy(filepath::AbstractString;
                            actor_lr::Float64=0.05, critic_lr::Float64=0.1,
                            gamma::Float64=0.0, temperature::Float64=1.0)
    open(filepath, "r") do f
        feat_dim = Int(read(f, Int32))
        num_actions = Int(read(f, Int32))
        Wpi = reshape(reinterpret(Float64, read(f, sizeof(Float64) * feat_dim * num_actions)),
                      feat_dim, num_actions)
        bpi = reinterpret(Float64, read(f, sizeof(Float64) * num_actions))
        Wv = reinterpret(Float64, read(f, sizeof(Float64) * feat_dim))
        bv = read(f, Float64)
        return GnnACPolicy(copy(Wpi), copy(bpi), copy(Wv), bv,
                           actor_lr, critic_lr, gamma, temperature, feat_dim, num_actions)
    end
end

function _softmax_logits(z::AbstractVector{Float64}, temperature::Float64)
    T = max(1e-6, temperature)
    zT = z ./ T
    m = maximum(zT)
    ex = exp.(zT .- m)
    s = sum(ex)
    return (s == 0.0) ? fill(1.0/length(z), length(z)) : (ex ./ s)
end

function _actor_logits(pol::GnnACPolicy, x::AbstractVector{Float64})
    return pol.Wpi' * x .+ pol.bpi
end

"""
    select_action_gnn(policy, x) -> Int (sampled)
"""
function select_action_gnn(policy::GnnACPolicy, x::AbstractVector{Float64})::Int
    probs = _softmax_logits(_actor_logits(policy, x), policy.temperature)
    r = rand()
    c = 0.0
    for i in 1:length(probs)
        c += probs[i]
        if r <= c
            return i
        end
    end
    return length(probs)
end

"""
    best_action_gnn(policy, x) -> Int (greedy)
"""
function best_action_gnn(policy::GnnACPolicy, x::AbstractVector{Float64})::Int
    probs = _softmax_logits(_actor_logits(policy, x), policy.temperature)
    return argmax(probs)
end

"""
    evaluate_gnn_actor_critic(env, policy, feature_fn) -> Float64
"""
function evaluate_gnn_actor_critic(env::CompressionEnv{T}, policy::GnnACPolicy,
                                   feature_fn) where {T<:Unsigned}
    state = reset!(env)
    while !env.done
        v = current_vertex(env)
        @assert v !== nothing
        x = feature_fn(v::T)
        a = best_action_gnn(policy, x)
        action = action_from_index(a)
        state, _, _ = step!(env, action)
    end
    return get_bits_per_edge(env)
end

"""
    train_gnn_actor_critic!(env, policy, feature_fn; episodes=50, eval_every=5, verbose=true)

Actor-critic with linear heads over GNN features.
feature_fn(v)::Vector{Float64} must return a feature vector of length feat_dim.
"""
function train_gnn_actor_critic!(env::CompressionEnv{T}, policy::GnnACPolicy, feature_fn;
                                 episodes::Int=50, eval_every::Int=5, verbose::Bool=true) where {T<:Unsigned}
    episode_bpe = Float64[]
    for ep in 1:episodes
        state = reset!(env)
        while !env.done
            v = current_vertex(env)
            @assert v !== nothing
            x = feature_fn(v::T)

            # Actor: pick action
            logits = _actor_logits(policy, x)
            probs = _softmax_logits(logits, policy.temperature)
            r = rand(); c = 0.0; a = 1
            for i in 1:length(probs)
                c += probs[i]
                if r <= c
                    a = i; break
                end
            end
            act = action_from_index(a)

            # Step environment
            next_state, reward, done = step!(env, act)

            # Critic: TD(0) baseline on x
            V_s = dot(policy.Wv, x) + policy.bv
            if done
                delta = reward - V_s
            else
                vnext = current_vertex(env)
                xnext = feature_fn(vnext::T)
                V_ns = dot(policy.Wv, xnext) + policy.bv
                delta = reward + policy.gamma * V_ns - V_s
            end
            # Update critic
            policy.Wv .+= policy.critic_lr * delta .* x
            policy.bv += policy.critic_lr * delta

            # Update actor: grad log π(a|x) = onehot(a) - probs
            for i in 1:policy.num_actions
                grad = ((i == a) ? 1.0 : 0.0) - probs[i]
                policy.Wpi[:, i] .+= policy.actor_lr * delta .* x .* grad
                policy.bpi[i] += policy.actor_lr * delta * grad
            end

            state = next_state
        end

        bpe = get_bits_per_edge(env)
        push!(episode_bpe, bpe)
        if verbose && (ep % eval_every == 0 || ep == 1)
            eval_bpe = evaluate_gnn_actor_critic(env, policy, feature_fn)
            println("  Episode $ep: train=$(round(bpe, digits=4)), eval=$(round(eval_bpe, digits=4))")
        end
    end
    return episode_bpe
end

"""
    train_gnn_ac_e2e!(env, policy, gnn, feature_fn; episodes=50, gnn_lr=0.0001,
                      gnn_update_every=5, eval_every=5, verbose=true)

End-to-end actor-critic training that fine-tunes GNN weights alongside linear heads.
Accumulates RL gradients w.r.t. GNN outputs (H1, Z2) over `gnn_update_every` episodes,
then backpropagates through the GNN and refreshes caches.
"""
function train_gnn_ac_e2e!(env::CompressionEnv{T}, policy::GnnACPolicy, gnn,
                           feature_fn;
                           episodes::Int=50, gnn_lr::Float64=0.0001,
                           gnn_update_every::Int=5, eval_every::Int=5,
                           verbose::Bool=true) where {T<:Unsigned}
    episode_bpe = Float64[]
    hidden_dim = size(gnn.H1, 2)

    # Accumulated gradients for GNN
    dL_dH1_acc = zeros(Float64, gnn.n, hidden_dim)
    dL_dZ2_acc = zeros(Float64, gnn.n)

    for ep in 1:episodes
        state = reset!(env)
        while !env.done
            v = current_vertex(env)
            @assert v !== nothing
            vi = Int(v)
            x = feature_fn(v::T)

            # Actor: pick action
            logits = _actor_logits(policy, x)
            probs = _softmax_logits(logits, policy.temperature)
            r = rand(); c = 0.0; a = 1
            for i in 1:length(probs)
                c += probs[i]
                if r <= c
                    a = i; break
                end
            end
            act = action_from_index(a)

            # Step environment
            next_state, reward, done = step!(env, act)

            # Critic: TD(0) baseline
            V_s = dot(policy.Wv, x) + policy.bv
            if done
                delta = reward - V_s
            else
                vnext = current_vertex(env)
                xnext = feature_fn(vnext::T)
                V_ns = dot(policy.Wv, xnext) + policy.bv
                delta = reward + policy.gamma * V_ns - V_s
            end

            # Update linear heads (same as train_gnn_actor_critic!)
            policy.Wv .+= policy.critic_lr * delta .* x
            policy.bv += policy.critic_lr * delta
            for i in 1:policy.num_actions
                grad = ((i == a) ? 1.0 : 0.0) - probs[i]
                policy.Wpi[:, i] .+= policy.actor_lr * delta .* x .* grad
                policy.bpi[i] += policy.actor_lr * delta * grad
            end

            # Accumulate gradient w.r.t. GNN features for this vertex
            # Actor gradient: dL/dx = δ · Σ_i (1{i==a} - π_i) · Wpi[:,i]
            dL_dx = zeros(Float64, policy.feat_dim)
            for i in 1:policy.num_actions
                grad_i = ((i == a) ? 1.0 : 0.0) - probs[i]
                dL_dx .+= delta * grad_i .* policy.Wpi[:, i]
            end
            # Critic gradient: dL/dx += δ · Wv
            dL_dx .+= delta .* policy.Wv

            # Split into H1 and Z2 components
            dL_dH1_acc[vi, :] .+= dL_dx[1:hidden_dim]
            dL_dZ2_acc[vi] += dL_dx[hidden_dim + 1]

            state = next_state
        end

        bpe = get_bits_per_edge(env)
        push!(episode_bpe, bpe)
        if verbose && (ep % eval_every == 0 || ep == 1)
            eval_bpe = evaluate_gnn_actor_critic(env, policy, feature_fn)
            println("  E2E Episode $ep: train=$(round(bpe, digits=4)), eval=$(round(eval_bpe, digits=4))")
        end

        # Periodically update GNN weights
        if ep % gnn_update_every == 0
            gnn_backward_rl!(gnn, dL_dZ2_acc, dL_dH1_acc; lr=gnn_lr)
            gnn_forward(gnn)  # refresh H1/Z2 caches
            fill!(dL_dH1_acc, 0.0)
            fill!(dL_dZ2_acc, 0.0)
            if verbose
                println("  [GNN update at episode $ep]")
            end
        end
    end
    return episode_bpe
end
