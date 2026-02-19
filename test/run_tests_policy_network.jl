include("run_tests_main.jl")

using .RL: PolicyNetwork, ActorCriticPN, train_ac_pn!, CompressionEnv, NUM_ACTIONS
using .GNN: GNNModel
using LightGraphs

@testset "Policy Network" begin
    # Create a simple graph
    g = SimpleDiGraph{UInt32}(4)
    add_edge!(g, 1, 2)
    add_edge!(g, 1, 3)
    add_edge!(g, 2, 4)
    add_edge!(g, 3, 4)

    # Create a compression environment
    env = CompressionEnv(g)

    # Create an actor-critic policy network
    ac = ActorCriticPN(5, 16, NUM_ACTIONS)

    # Capture initial weights
    initial_weights = copy(ac.policy_network.model[1].weight)

    # Train the policy network
    train_ac_pn!(env, ac, 10)

    # Test that the policy network weights have changed
    @test ac.policy_network.model[1].weight != initial_weights
    @test ac.policy_network isa PolicyNetwork
end

@testset "GAT Network" begin
    # Create a simple graph
    g = SimpleDiGraph{UInt32}(4)
    add_edge!(g, 1, 2)
    add_edge!(g, 1, 3)
    add_edge!(g, 2, 4)
    add_edge!(g, 3, 4)

    # Create a GNN model with GAT enabled
    gnn = GNNModel(g, use_gat=true)
    
    # Check that the GAT layers were created
    @test gnn.use_gat == true
    @test length(gnn.gat_layers) == 2

    # More tests to come here once the GAT backward pass is fully implemented
end
