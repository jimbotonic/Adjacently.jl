"""
    Agent

Per-agent state for the Mycelial Polis model (roadmap §6).

Fields map to the §6 state vector `x_i(t) = (a_i, c_i, r_i, l_i, s_i, q_i,
e_i, b_i, z_i)` plus the §11 split of legibility into internal and external
components. Trust/reputation `q_i` and exposure `e_i` are derived on demand
from the multiplex layers and not stored here.

Roles (stage `z_i`, §6, §8): `:outsider, :observer, :sympathizer, :user,
:contributor, :steward, :removed`. Adoption advances by stage thresholds in
`dynamics.jl`. The `:removed` role represents a host node-removal attack
(arrest, ban, account deletion).
"""
Base.@kwdef mutable struct Agent
    id::Int
    awareness::Float32   = 0f0   # a_i — adoption level in [0,1]
    commitment::Float32  = 0f0   # c_i — durability of participation
    capacity::Float32    = 1f0   # s_i — resource / contribution capacity
    fear::Float32        = 0f0   # r_i — perceived repression risk
    backfire::Float32    = 0f0   # b_i — outrage / solidarity boost
    identity::Float32    = 0f0   # χ·I_i — identity alignment with the polis
    L_int::Float32       = 0f0   # internal legibility (§11)
    L_ext::Float32       = 0f0   # external legibility (§11)
    role::Symbol         = :outsider
    # Paper 2 — E2 endogenous conflict.
    # `faction`: which internal coalition the agent belongs to;
    # `:none` until a disagreement event labels them. Capped at 4
    # distinct factions to keep state bounded.
    # `disagreement`: per-agent accumulator that triggers a schism
    # (faction flip + trust-edge decay) when it crosses
    # `params.schism_threshold`. Decays at `params.disagreement_decay`.
    faction::Symbol      = :none
    disagreement::Float32 = 0f0
end

const ROLE_RANK = Dict{Symbol,Int}(
    :outsider     => 0,
    :observer     => 1,
    :sympathizer  => 2,
    :user         => 3,
    :contributor  => 4,
    :steward      => 5,
    :removed      => -1,
    :infiltrator  => -2,    # in the network, not on the polis side (item 5)
    :defector     => -3,    # was committed; abandoned (item 6 — wires triggers)
)

is_active(a::Agent)      = a.role !== :removed && a.role !== :outsider
is_committed(a::Agent)   = ROLE_RANK[a.role] >= ROLE_RANK[:user]
is_steward(a::Agent)     = a.role === :steward
is_infiltrator(a::Agent) = a.role === :infiltrator
is_defector(a::Agent)    = a.role === :defector
