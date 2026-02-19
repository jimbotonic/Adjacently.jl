#
# GNN-augmented actor-critic: linear heads over per-vertex GNN features.
#

using LinearAlgebra: dot, norm
using BSON: @save, @load
using Statistics: mean

mutable struct GnnACPolicy
    Wpi::Matrix{Float64}
    bpi::Vector{Float64}
    Wv::Vector{Float64}
    bv::Float64
    actor_lr::Float64
    critic_lr::Float64
    gamma::Float64
    temperature::Float64
    entropy_coeff::Float64
    feat_dim::Int
    num_actions::Int
end

function GnnACPolicy(feat_dim::Int; num_actions::Int=NUM_ACTIONS,
                     actor_lr::Float64=0.01, critic_lr::Float64=0.01,
                     gamma::Float64=0.0, temperature::Float64=2.0,
                     entropy_coeff::Float64=0.01)
    Wpi = (rand(Float64, feat_dim, num_actions) .- 0.5) .* 0.01
    bpi = zeros(Float64, num_actions)
    Wv = (rand(Float64, feat_dim) .- 0.5) .* 0.01
    bv = 0.0
    return GnnACPolicy(Wpi, bpi, Wv, bv, actor_lr, critic_lr, gamma, temperature, entropy_coeff, feat_dim, num_actions)
end

function _softmax_logits(z::AbstractVector{Float64}, temperature::Float64)
    T = max(1e-6, temperature)
    zT = z ./ T
    m = maximum(zT)
    ex = exp.(zT .- m)
    s = sum(ex)
    if s <= 0.0 || isnan(s)
        return fill(1.0/length(z), length(z))
    end
    return ex ./ s
end

function _actor_logits(pol::GnnACPolicy, x::AbstractVector{Float64})
    logits = pol.Wpi' * x .+ pol.bpi
    clamp!(logits, -20.0, 20.0)
    return logits
end

function select_action_gnn(policy::GnnACPolicy, x::AbstractVector{Float64})::Int
    logits = _actor_logits(policy, x)
    probs = _softmax_logits(logits, policy.temperature)
    r = rand(); c = 0.0
    for i in 1:length(probs)
        c += probs[i]
        if r <= c; return i; end
    end
    return length(probs)
end

function best_action_gnn(policy::GnnACPolicy, x::AbstractVector{Float64})::Int
    logits = _actor_logits(policy, x)
    probs = _softmax_logits(logits, policy.temperature)
    return argmax(probs)
end

function evaluate_gnn_actor_critic(env, policy::GnnACPolicy, feature_fn)
    reset!(env)
    while !env.done
        v = current_vertex(env)
        x = feature_fn(v)
        a = best_action_gnn(policy, x)
        step!(env, action_from_index(a))
    end
    return get_bits_per_edge(env)
end

function train_gnn_ac_e2e!(env, policy::GnnACPolicy, gnn,
                           feature_fn;
                           episodes::Int=50, gnn_lr::Float64=0.0001,
                           gnn_update_every::Int=5, eval_every::Int=5,
                           verbose::Bool=true)
    episode_bpe = Float64[]
    hidden_dim = size(gnn.H1, 2)
    dL_dH1_acc = zeros(Float64, gnn.n, hidden_dim)
    dL_dZ2_acc = zeros(Float64, gnn.n)
    action_counts = zeros(Int, policy.num_actions)

    for ep in 1:episodes
        state = reset!(env)
        fill!(action_counts, 0)
        ep_reward = 0.0; n_steps = 0
        deltas = Float64[]
        entropies = Float64[]
        
        while !env.done
            v = current_vertex(env)
            @assert v !== nothing
            vi = Int(v); x = feature_fn(v)
            if any(isnan, x); x = replace(x, NaN => 0.0); end

            logits = _actor_logits(policy, x)
            probs = _softmax_logits(logits, policy.temperature)
            log_probs = log.(probs .+ 1e-9)
            
            # Entropy for logging
            ent = -sum(probs .* log_probs)
            push!(entropies, ent)
            
            r = rand(); c = 0.0; a = 1
            for i in 1:length(probs)
                c += probs[i]
                if r <= c; a = i; break; end
            end
            action_counts[a] += 1
            act = action_from_index(a)

            next_state, reward, done = step!(env, act)
            ep_reward += reward; n_steps += 1

            V_s = dot(policy.Wv, x) + policy.bv
            if done
                delta = reward - V_s
            else
                vnext = current_vertex(env)
                xnext = feature_fn(vnext)
                V_ns = dot(policy.Wv, xnext) + policy.bv
                delta = reward + policy.gamma * V_ns - V_s
            end
            
            if isnan(delta); delta = 0.0; end
            delta = clamp(delta, -50.0, 50.0)
            push!(deltas, delta)

            # Update heads
            policy.Wv .+= policy.critic_lr * delta .* x
            policy.bv += policy.critic_lr * delta

            for i in 1:policy.num_actions
                pg_grad = ((i == a) ? 1.0 : 0.0) - probs[i]
                # Entropy gradient: push logits away from being peaked
                ent_grad = -(log_probs[i] - sum(probs .* log_probs)) * probs[i]
                combined_grad = delta * pg_grad + policy.entropy_coeff * ent_grad
                if !isnan(combined_grad)
                    policy.Wpi[:, i] .+= policy.actor_lr .* x .* combined_grad
                    policy.bpi[i] += policy.actor_lr * combined_grad
                end
            end

            dL_dx = zeros(Float64, policy.feat_dim)
            for i in 1:policy.num_actions
                grad_i = ((i == a) ? 1.0 : 0.0) - probs[i]
                dL_dx .+= delta * grad_i .* policy.Wpi[:, i]
            end
            dL_dx .+= delta .* policy.Wv
            
            if !any(isnan, dL_dx)
                dL_dH1_acc[vi, :] .+= dL_dx[1:hidden_dim]
                dL_dZ2_acc[vi] += dL_dx[hidden_dim + 1]
            end
            state = next_state
        end

        bpe = get_bits_per_edge(env)
        push!(episode_bpe, bpe)
        if verbose && (ep % eval_every == 0 || ep == 1)
            eval_bpe = evaluate_gnn_actor_critic(env, policy, feature_fn)
            avg_rew = ep_reward / n_steps
            avg_delta = mean(abs.(deltas))
            avg_ent = mean(entropies)
            most_common_a = argmax(action_counts)
            most_common_pct = 100.0 * action_counts[most_common_a] / n_steps
            println("  Ep $ep: train=$(round(bpe, digits=4)), eval=$(round(eval_bpe, digits=4)), rew=$(round(avg_rew, digits=2)), delta=$(round(avg_delta, digits=2)), ent=$(round(avg_ent, digits=2)), top_act=$most_common_a ($(round(most_common_pct, digits=1))%)")
        end

        if ep % gnn_update_every == 0
            n_denom = max(gnn.n * gnn_update_every, 1)
            dL_dH1_acc ./= n_denom; dL_dZ2_acc ./= n_denom
            gnn_grad_norm = norm(dL_dH1_acc) + norm(dL_dZ2_acc)
            if !isnan(gnn_grad_norm)
                if gnn_grad_norm > 1.0
                    scale = 1.0 / gnn_grad_norm
                    dL_dH1_acc .*= scale; dL_dZ2_acc .*= scale
                end
                gnn_backward_rl!(gnn, dL_dZ2_acc, dL_dH1_acc; lr=gnn_lr)
                gnn_forward(gnn)
            end
            fill!(dL_dH1_acc, 0.0); fill!(dL_dZ2_acc, 0.0)
        end
    end
    return episode_bpe
end

function save_gnn_ac_policy(policy::GnnACPolicy, gnn::GNNModel, filepath::AbstractString)
    @save filepath policy gnn
end

function load_gnn_ac_policy(filepath::AbstractString;
                            actor_lr::Float64=0.05, critic_lr::Float64=0.1,
                            gamma::Float64=0.0, temperature::Float64=1.0)
    @load filepath policy gnn
    policy.actor_lr = actor_lr; policy.critic_lr = critic_lr
    policy.gamma = gamma; policy.temperature = temperature
    return policy, gnn
end
