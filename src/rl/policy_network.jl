#
# Policy Network for RL-based compression policy learning.
#

using Flux
using BSON: @save, @load
using StatsBase
using Optimisers

mutable struct PolicyNetwork
    model::Chain
    opt_state::Any
end

function PolicyNetwork(input_dim::Int, hidden_dim::Int, num_actions::Int; lr=0.001)
    model = Chain(
        Dense(input_dim, hidden_dim, relu),
        Dense(hidden_dim, num_actions),
        softmax
    )
    opt_state = Optimisers.setup(Optimisers.Adam(lr), model)
    return PolicyNetwork(model, opt_state)
end

function select_action_pn(pn::PolicyNetwork, state)
    probs = pn.model(state)
    a = sample(1:length(probs), Weights(probs))
    return a, log(probs[a] + 1f-8)
end

function train_pn!(pn::PolicyNetwork, state, action, advantage)
    function loss(m)
        probs = m(state)
        lp = log(probs[action] + 1f-8)
        return -lp * advantage
    end
    grads = Flux.gradient(loss, pn.model)[1]
    pn.opt_state, pn.model = Optimisers.update(pn.opt_state, pn.model, grads)
end

function save_pn_policy(pn::PolicyNetwork, filepath::String)
    @save filepath pn
end

function load_pn_policy(filepath::String)::PolicyNetwork
    @load filepath pn
    return pn
end
