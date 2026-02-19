#
# Tabular Q-learning policy for compression encoding decisions.
# Uses epsilon-greedy exploration with decay.
#

mutable struct QPolicy
    Q::Matrix{Float64}        # (num_states, num_actions)
    alpha::Float64             # learning rate
    gamma::Float64             # discount factor
    epsilon::Float64           # exploration rate
    epsilon_decay::Float64     # decay per episode
    epsilon_min::Float64       # minimum epsilon
    num_states::Int
    num_actions::Int
    action_counts::Matrix{Int} # track action usage per state
end

function QPolicy(; num_states::Int=NUM_STATES,
                   num_actions::Int=NUM_ACTIONS,
                   alpha::Float64=0.1,
                   gamma::Float64=0.99,
                   epsilon::Float64=0.3,
                   epsilon_decay::Float64=0.95,
                   epsilon_min::Float64=0.05)
    Q = zeros(Float64, num_states, num_actions)
    action_counts = zeros(Int, num_states, num_actions)
    return QPolicy(Q, alpha, gamma, epsilon, epsilon_decay, epsilon_min,
                   num_states, num_actions, action_counts)
end

"""
    select_action(policy, state_idx) -> Int

Epsilon-greedy action selection. Returns action index (1-based).
"""
function select_action(policy::QPolicy, state_idx::Int)::Int
    if rand() < policy.epsilon
        return rand(1:policy.num_actions)
    else
        return best_action(policy, state_idx)
    end
end

"""
    best_action(policy, state_idx) -> Int

Return the action with highest Q-value for the given state (greedy).
"""
function best_action(policy::QPolicy, state_idx::Int)::Int
    row = view(policy.Q, state_idx, :)
    return argmax(row)
end

"""
    update!(policy, state, action, reward, next_state)

Q-learning update: Q(s,a) += alpha * (r + gamma * max_a' Q(s',a') - Q(s,a))
"""
function update!(policy::QPolicy, state::Int, action::Int,
                 reward::Float64, next_state::Int)
    current_q = policy.Q[state, action]
    max_next_q = maximum(view(policy.Q, next_state, :))
    td_target = reward + policy.gamma * max_next_q
    policy.Q[state, action] = current_q + policy.alpha * (td_target - current_q)
    policy.action_counts[state, action] += 1
end

"""
    decay_epsilon!(policy)

Decay exploration rate after each episode.
"""
function decay_epsilon!(policy::QPolicy)
    policy.epsilon = max(policy.epsilon_min, policy.epsilon * policy.epsilon_decay)
end

"""
    get_policy_summary(policy) -> String

Return a summary of the most-chosen actions per state bucket.
"""
function get_policy_summary(policy::QPolicy)::String
    lines = String[]
    push!(lines, "Policy Summary (top actions per state category):")
    push!(lines, "-" ^ 60)

    # Summarize by degree bin (aggregate across other features)
    degree_labels = ["deg=0", "deg=1-3", "deg=4-10", "deg=11-50", "deg=51-200", "deg=201+"]
    for d in 1:DEGREE_BINS
        # Find all states with this degree bin
        total_counts = zeros(Int, policy.num_actions)
        for iv in 1:INTERVAL_DENSITY_BINS
            for g in 1:MAX_GAP_BINS
                for o in 1:REF_OVERLAP_BINS
                    for rw in 1:REF_WINDOW_DENSITY_BINS
                        f = VertexFeatures(d, iv, g, o, rw)
                        s = feature_index(f)
                        if 1 <= s <= policy.num_states
                            total_counts .+= policy.action_counts[s, :]
                        end
                    end
                end
            end
        end

        total_visits = sum(total_counts)
        total_visits == 0 && continue

        # Top 3 actions
        top_idxs = sortperm(total_counts, rev=true)[1:min(3, policy.num_actions)]
        top_strs = String[]
        for idx in top_idxs
            total_counts[idx] == 0 && continue
            a = action_from_index(idx)
            pct = round(100.0 * total_counts[idx] / total_visits, digits=1)
            push!(top_strs, "$(a.reference_mode)/$(a.encoding_type)/mil=$(a.min_interval_length) ($(pct)%)")
        end

        push!(lines, "  $(degree_labels[d]) [$(total_visits) visits]: $(join(top_strs, ", "))")
    end

    return join(lines, "\n")
end

const QPOLICY_MAGIC = UInt8[0x52, 0x4c, 0x51, 0x50]  # "RLQP"
const QPOLICY_VERSION = UInt16(1)

"""
    save_policy(policy::QPolicy, filepath::String)

Save a trained Q-policy to a binary file (.qpolicy format).

Format:
- Header (16 bytes): magic "RLQP" (4) + version UInt16 (2) + num_states UInt32 (4) + num_actions UInt32 (4) + reserved UInt16 (2)
- Hyperparameters (40 bytes): alpha, gamma, epsilon, epsilon_decay, epsilon_min as Float64
- Q-matrix (num_states * num_actions * 8 bytes): Float64 row-major
- Action counts (num_states * num_actions * 8 bytes): Int64 row-major
"""
function save_policy(policy::QPolicy, filepath::String)
    open(filepath, "w") do f
        # Header
        write(f, QPOLICY_MAGIC)
        write(f, QPOLICY_VERSION)
        write(f, UInt32(policy.num_states))
        write(f, UInt32(policy.num_actions))
        write(f, UInt16(0))  # reserved

        # Hyperparameters
        write(f, policy.alpha)
        write(f, policy.gamma)
        write(f, policy.epsilon)
        write(f, policy.epsilon_decay)
        write(f, policy.epsilon_min)

        # Q-matrix (row-major)
        for s in 1:policy.num_states
            for a in 1:policy.num_actions
                write(f, policy.Q[s, a])
            end
        end

        # Action counts (row-major)
        for s in 1:policy.num_states
            for a in 1:policy.num_actions
                write(f, Int64(policy.action_counts[s, a]))
            end
        end
    end
    return nothing
end

"""
    load_policy(filepath::String) -> QPolicy

Load a Q-policy from a binary .qpolicy file.
"""
function load_policy(filepath::String)::QPolicy
    open(filepath, "r") do f
        # Header
        magic = read(f, 4)
        if magic != QPOLICY_MAGIC
            error("Invalid policy file: bad magic (expected RLQP)")
        end
        version = read(f, UInt16)
        if version != QPOLICY_VERSION
            error("Unsupported policy version: $version (expected $QPOLICY_VERSION)")
        end
        num_states = Int(read(f, UInt32))
        num_actions = Int(read(f, UInt32))
        _ = read(f, UInt16)  # reserved

        # Hyperparameters
        alpha = read(f, Float64)
        gamma = read(f, Float64)
        epsilon = read(f, Float64)
        epsilon_decay = read(f, Float64)
        epsilon_min = read(f, Float64)

        # Q-matrix
        Q = zeros(Float64, num_states, num_actions)
        for s in 1:num_states
            for a in 1:num_actions
                Q[s, a] = read(f, Float64)
            end
        end

        # Action counts
        action_counts = zeros(Int, num_states, num_actions)
        for s in 1:num_states
            for a in 1:num_actions
                action_counts[s, a] = Int(read(f, Int64))
            end
        end

        return QPolicy(Q, alpha, gamma, epsilon, epsilon_decay, epsilon_min,
                       num_states, num_actions, action_counts)
    end
end
