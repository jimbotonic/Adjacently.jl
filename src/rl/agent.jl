#
# Training loop and evaluation for RL-based compression policy learning.
#

struct TrainingConfig
    num_episodes::Int
    eval_every::Int
    verbose::Bool
end

TrainingConfig(; num_episodes::Int=50, eval_every::Int=5, verbose::Bool=true) =
    TrainingConfig(num_episodes, eval_every, verbose)

struct TrainingResult
    episode_bpe::Vector{Float64}
    best_bpe::Float64
    best_episode::Int
    baseline_bpe::Float64
    policy::QPolicy
end

"""
    run_baseline(env) -> Float64

Run one episode with the default greedy action (fibonacci + reference + interval=4).
Returns bits per edge.
"""
function run_baseline(env::CompressionEnv{T})::Float64 where {T<:Unsigned}
    default_action = Action(:fibonacci, :reference, 4)
    state = reset!(env)
    while !env.done
        _, _, _ = step!(env, default_action)
    end
    return get_bits_per_edge(env)
end

"""
    evaluate(env, policy) -> Float64

Run one episode using the greedy policy (no exploration). Returns bits per edge.
"""
function evaluate(env::CompressionEnv{T}, policy::QPolicy)::Float64 where {T<:Unsigned}
    state = reset!(env)
    while !env.done
        s_idx = feature_index(state)
        a_idx = best_action(policy, s_idx)
        action = action_from_index(a_idx)
        state, _, done = step!(env, action)
    end
    return get_bits_per_edge(env)
end

"""
    train!(env, policy, config) -> TrainingResult

Main training loop. Runs Q-learning episodes on the compression environment.
"""
function train!(env::CompressionEnv{T}, policy::QPolicy,
                config::TrainingConfig=TrainingConfig()) where {T<:Unsigned}

    # Baseline measurement
    if config.verbose
        println("Running baseline (fibonacci + reference + interval=4)...")
    end
    baseline_bpe = run_baseline(env)
    if config.verbose
        println("  Baseline: $(round(baseline_bpe, digits=4)) bits/edge")
        println("")
        println("Starting RL training ($(config.num_episodes) episodes)...")
        println("-" ^ 60)
    end

    episode_bpe = Float64[]
    best_bpe = Inf
    best_episode = 0

    for ep in 1:config.num_episodes
        state = reset!(env)

        while !env.done
            s_idx = feature_index(state)
            a_idx = select_action(policy, s_idx)
            action = action_from_index(a_idx)

            next_state, reward, done = step!(env, action)

            # Q-learning update
            ns_idx = done ? s_idx : feature_index(next_state)
            update!(policy, s_idx, a_idx, reward, ns_idx)

            state = next_state
        end

        decay_epsilon!(policy)
        bpe = get_bits_per_edge(env)
        push!(episode_bpe, bpe)

        if bpe < best_bpe
            best_bpe = bpe
            best_episode = ep
        end

        if config.verbose && (ep % config.eval_every == 0 || ep == 1)
            eval_bpe = evaluate(env, policy)
            println("  Episode $ep: train=$(round(bpe, digits=4)), eval=$(round(eval_bpe, digits=4)), " *
                    "best=$(round(best_bpe, digits=4)), eps=$(round(policy.epsilon, digits=3))")
        end
    end

    # Final evaluation
    if config.verbose
        println("-" ^ 60)
        final_bpe = evaluate(env, policy)
        improvement = baseline_bpe - final_bpe
        pct = 100.0 * improvement / baseline_bpe
        println("")
        println("Training complete!")
        println("  Baseline:       $(round(baseline_bpe, digits=4)) bits/edge")
        println("  Best (train):   $(round(best_bpe, digits=4)) bits/edge (episode $best_episode)")
        println("  Final (greedy): $(round(final_bpe, digits=4)) bits/edge")
        println("  Improvement:    $(round(improvement, digits=4)) bits/edge ($(round(pct, digits=2))%)")
        println("")
        println(get_policy_summary(policy))
    end

    return TrainingResult(episode_bpe, best_bpe, best_episode, baseline_bpe, policy)
end
