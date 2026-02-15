#
# On-policy actor-critic (tabular softmax policy + state-value baseline)
# Optimizes per-vertex encoding actions to reduce bits/edge.
#

mutable struct ACPolicy
    logits::Matrix{Float64}   # (num_states, num_actions): unnormalized preferences
    V::Vector{Float64}        # (num_states): state-value baseline
    actor_lr::Float64
    critic_lr::Float64
    gamma::Float64
    temperature::Float64
    num_states::Int
    num_actions::Int
end

function ACPolicy(; num_states::Int=NUM_STATES,
                    num_actions::Int=NUM_ACTIONS,
                    actor_lr::Float64=0.05,
                    critic_lr::Float64=0.1,
                    gamma::Float64=0.0,           # episodic, per-step rewards sum; gamma=0 is fine
                    temperature::Float64=1.0)
    logits = zeros(Float64, num_states, num_actions)
    V = zeros(Float64, num_states)
    return ACPolicy(logits, V, actor_lr, critic_lr, gamma, temperature,
                    num_states, num_actions)
end

"""
    _softmax_row(logits_row, temperature) -> prob vector
"""
function _softmax_row(row::AbstractVector{Float64}, temperature::Float64)
    T = max(1e-6, temperature)
    z = row ./ T
    m = maximum(z)
    ex = exp.(z .- m)
    s = sum(ex)
    return (s == 0.0) ? fill(1.0/length(row), length(row)) : (ex ./ s)
end

"""
    select_action_ac(policy, state_idx) -> Int (sampled)
"""
function select_action_ac(policy::ACPolicy, state_idx::Int)::Int
    probs = _softmax_row(view(policy.logits, state_idx, :), policy.temperature)
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
    best_action_ac(policy, state_idx) -> Int (greedy)
"""
function best_action_ac(policy::ACPolicy, state_idx::Int)::Int
    probs = _softmax_row(view(policy.logits, state_idx, :), policy.temperature)
    return argmax(probs)
end

"""
    train_actor_critic!(env, policy; episodes=50, eval_every=5, verbose=true)

On-policy actor-critic with TD(0) baseline and softmax policy.
Returns vector of episode bpe measurements.
"""
function train_actor_critic!(env::CompressionEnv{T}, policy::ACPolicy;
                             episodes::Int=50, eval_every::Int=5, verbose::Bool=true) where {T<:Unsigned}
    episode_bpe = Float64[]

    if verbose
        println("Running actor-critic training ($episodes episodes)...")
        println("-" ^ 60)
    end

    for ep in 1:episodes
        state = reset!(env)
        while !env.done
            s = feature_index(state)
            a = select_action_ac(policy, s)
            act = action_from_index(a)
            next_state, reward, done = step!(env, act)

            # TD(0) update for critic
            ns = done ? s : feature_index(next_state)
            td_target = reward + policy.gamma * policy.V[ns]
            delta = td_target - policy.V[s]
            policy.V[s] += policy.critic_lr * delta

            # Policy gradient step: grad log pi(s,a) * advantage (delta)
            probs = _softmax_row(view(policy.logits, s, :), policy.temperature)
            for i in 1:policy.num_actions
                grad = (i == a ? 1.0 : 0.0) - probs[i]
                policy.logits[s, i] += policy.actor_lr * delta * grad
            end

            state = next_state
        end

        bpe = get_bits_per_edge(env)
        push!(episode_bpe, bpe)

        if verbose && (ep % eval_every == 0 || ep == 1)
            eval_bpe = evaluate_actor_critic(env, policy)
            println("  Episode $ep: train=$(round(bpe, digits=4)), eval=$(round(eval_bpe, digits=4))")
        end
    end

    return episode_bpe
end

"""
    evaluate_actor_critic(env, policy) -> Float64

Evaluate the greedy actor policy (argmax over softmax probs).
"""
function evaluate_actor_critic(env::CompressionEnv{T}, policy::ACPolicy)::Float64 where {T<:Unsigned}
    state = reset!(env)
    while !env.done
        s = feature_index(state)
        a = best_action_ac(policy, s)
        action = action_from_index(a)
        state, _, _ = step!(env, action)
    end
    return get_bits_per_edge(env)
end

