# Bounded, deterministic search/emission phases. Search owns no output stream.
# Each worker owns its workspace, even if a task migrates between Julia threads.
function _check_parallel_search(enabled, coding_scheme, exact_costing, workers, batch_size)
    enabled || return nothing
    coding_scheme == :children || throw(ArgumentError("parallel_search currently supports children mode only"))
    exact_costing && throw(ArgumentError("parallel_search supports FULL/FAST analytical cost models, not exact_costing trial writes"))
    workers > 0 || throw(ArgumentError("search_workers must be positive"))
    batch_size > 0 || throw(ArgumentError("search_batch_size must be positive"))
    return nothing
end

function _parallel_vertex_batches!(search::F, emit::E,
        neighbor_lists::Dict{T,Vector{T}}, vs::Int, ref_window_size::Int;
        workers::Int, batch_size::Int) where {T<:Unsigned,F,E}
    ref_window_size >= 0 || throw(ArgumentError("reference window must be nonnegative"))
    vs == 0 && return nothing
    nworkers = min(workers, Threads.nthreads(:default), batch_size, vs)
    states = [(VertexSearchWorkspace{T}(), T[]) for _ in 1:nworkers]
    for first in 1:batch_size:vs
        count = min(batch_size, vs-first+1)
        lists = Vector{Vector{T}}(undef, count)
        choices = Vector{Any}(undef, count)
        next = Threads.Atomic{Int}(1)
        function search_batch!(state)
            ws, window = state
            while true
                i = Threads.atomic_add!(next, 1)
                i > count && break
                v = first+i-1
                neighbors = sort(get(neighbor_lists, T(v), T[]))
                lists[i] = neighbors
                # Serial search includes empty earlier lists in its window too.
                start = max(1, v-ref_window_size)
                resize!(window, v-start)
                for j in eachindex(window); window[j] = T(start+j-1); end
                choices[i] = isempty(neighbors) ? nothing : search(T(v), neighbors, window, ws)
            end
        end
        if nworkers == 1
            search_batch!(states[1])
        else
            @sync for state in states
                Threads.@spawn search_batch!(state)
            end
        end
        # Barrier above: no cost search overlaps serial adaptive coder/tap writes.
        # Search results own their winning bitmaps/residuals, not workspace views.
        for i in 1:count
            emit(T(first+i-1), lists[i], choices[i])
        end
    end
    return nothing
end
