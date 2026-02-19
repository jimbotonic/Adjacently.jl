#
# Actor-Critic learning with a policy network.
#

using Flux
using StatsBase
using Optimisers

mutable struct ActorCriticPN
    policy_network::PolicyNetwork
    value_network::Chain
    opt_state::Any
    gamma::Float64
end

function ActorCriticPN(input_dim::Int, hidden_dim::Int, num_actions::Int; lr=0.001, gamma=0.99)
    policy_network = PolicyNetwork(input_dim, hidden_dim, num_actions; lr=lr)
    value_network = Chain(
        Dense(input_dim, hidden_dim, relu),
        Dense(hidden_dim, 1)
    )
    # Combined model for optimization
    model = (pn = policy_network.model, vn = value_network)
    opt_state = Optimisers.setup(Optimisers.Adam(lr), model)
    return ActorCriticPN(policy_network, value_network, opt_state, gamma)
end

function train_ac_pn!(env::CompressionEnv, ac::ActorCriticPN, episodes::Int)
    for ep in 1:episodes
        state = reset!(env)
        done = false
        while !done
            state_vec = Float32[
                state.degree_bin, 
                state.interval_density, 
                state.max_gap_bin, 
                state.ref_overlap_bin, 
                state.ref_window_density_bin
            ]
            
            action, _ = select_action_pn(ac.policy_network, state_vec)
            
            next_state, reward, done = step!(env, action_from_index(action))
            
            next_state_vec = Float32[
                next_state.degree_bin,
                next_state.interval_density,
                next_state.max_gap_bin,
                next_state.ref_overlap_bin,
                next_state.ref_window_density_bin
            ]
            
            # Loss function for explicit gradients
            function loss_fn(m)
                v_s = m.vn(state_vec)[1]
                # Standard AC target: reward + gamma * V(s')
                # We use ignore_derivatives for the next state value as is standard in TD
                v_ns = Flux.ignore_derivatives() do 
                    m.vn(next_state_vec)[1]
                end
                
                target = reward + ac.gamma * v_ns * (1 - done)
                advantage = target - v_s
                
                # Re-calculate log_prob for the action taken to allow backprop
                probs = m.pn(state_vec)
                lp = log(probs[action] + 1f-8)
                
                actor_loss = -lp * Flux.ignore_derivatives(() -> advantage)
                critic_loss = (target - v_s)^2
                
                return actor_loss + critic_loss
            end
            
            model = (pn = ac.policy_network.model, vn = ac.value_network)
            grads = Flux.gradient(loss_fn, model)[1]
            ac.opt_state, model = Optimisers.update(ac.opt_state, model, grads)
            
            # Update models back into ac
            ac.policy_network.model = model.pn
            ac.value_network = model.vn
            
            state = next_state
        end
    end
end
