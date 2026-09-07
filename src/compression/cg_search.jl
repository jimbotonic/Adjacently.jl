# CG cost search owns no stream. Worker scratch is task-owned, not thread-indexed.
function _check_parallel_search(params, enabled, workers)
    enabled || return nothing
    workers > 0 || throw(ArgumentError("search_workers must be positive"))
    (params.intra_mgs || params.intra_block_try) && throw(ArgumentError(
        "CG parallel_search supports the per-vertex analytical path, not MGS/block trial encoders"))
    return nothing
end

function _search_cluster(all_nl::Vector{Vector{T}}, params;
        parallel_search=false, search_workers=Threads.nthreads(:default),
        progress=nothing, ci=1, n_clusters=1) where {T<:Unsigned}
    s = length(all_nl)
        # The decision arrays were already cluster-sized in the serial encoder.
    use_ref_vec = Vector{Bool}(undef, s)
    ref_delta_vec = Vector{UInt32}(undef, s)
    ref_positions_list = Vector{Vector{T}}(undef, s)
    additions_list = Vector{Vector{T}}(undef, s)
    ref_len_list = Vector{Int}(undef, s)
    # Decisions for encoding modes (stored during search, replayed in write phase)
    copy_mode_vec = Vector{UInt8}(undef, s)    # 0=bitmap, 1=copy-blocks, 2=complement
    add_mode_vec = Vector{Bool}(undef, s)      # true=intervals, false=stop-delta (ref adds)
    raw_mode_vec = Vector{Bool}(undef, s)      # true=intervals, false=stop-delta (raw)
    # Per-vertex mil (greedy search populates this; otherwise fixed)
    mil_vec = fill(params.intra_mil, s)
    _mil_options = [2, 3, 4, 5]

    function search_vertices!(next)
        # Double-buffer swap: two position/adds vector pairs (optimization #5)
        _pos_a = T[]; _adds_a = T[]
        _pos_b = T[]; _adds_b = T[]

        idx_local = 0
        while true
            idx_local = next === nothing ? idx_local + 1 : Threads.atomic_add!(next, 1)
            idx_local > s && break
            next === nothing && progress !== nothing && progress(idx_local, s, ci, n_clusters)
            nl = all_nl[idx_local]
            # decide reference and mil
            use_ref = false; ref_delta_val = UInt32(0)
            best_copy_mode = 0x00; best_add_mode = false; best_raw_mode = false
            if params.intra_greedy_mil
                # Greedy per-vertex mil search: try all mil values for raw and ref (analytical)
                best_bits = typemax(Int)
                best_mil_val = params.intra_mil
                best_is_ref = false
                best_ref_idx = 0

                _zz_vid = params.intra_zigzag ? T(idx_local) : nothing
                if params.cost_model == Compression.COST_MODEL_FAST
                    # Fast model: single MIL, compare interval vs stop-delta
                    fast_mil = params.intra_adapt_mil > 0 ? params.intra_adapt_mil : params.intra_mil
                    iv_raw = estimate_interval_runlength_encoding_cost(nl, :fibonacci, fast_mil, 3; vertex_id=_zz_vid)
                    sd_raw = if params.intra_stop_deltas
                        _estimate_stop_delta_zigzag_cost(nl, :fibonacci, _zz_vid)
                    else
                        ab = _estimate_small_count_cost(length(nl), params.count_varint)
                        if !isempty(nl)
                            ab += _estimate_delta_list_cost(nl, :fibonacci; vertex_id=_zz_vid)
                        end
                        ab
                    end
                    raw_bits = min(iv_raw, sd_raw)
                    if raw_bits < best_bits
                        best_bits = raw_bits
                        best_mil_val = fast_mil
                        best_is_ref = false
                    end
                else
                for mil in _mil_options
                    if params.intra_lr_split
                        raw_bits = _estimate_ir_lr_cost(nl, :fibonacci, mil, _zz_vid; tight_deltas=params.intra_tight_deltas)
                    else
                        raw_bits = estimate_interval_runlength_encoding_cost(nl, :fibonacci, mil, 3; vertex_id=_zz_vid)
                    end
                    if raw_bits < best_bits
                        best_bits = raw_bits
                        best_mil_val = mil
                        best_is_ref = false
                    end
                end
                end

                # Try ref encoding with 2-phase pruning + analytical cost
                if params.intra_ref_enabled && idx_local > 1
                    wstart = max(1, idx_local - params.intra_ref_window)
                    wend = idx_local - 1
                    n_candidates = wend - wstart + 1

                    # Phase 1: overlap screening (cheap O(|nl|+|ref|) per candidate)
                    _max_k_greedy = params.cost_model == Compression.COST_MODEL_FAST ? MAX_REF_CANDIDATES_PHASE2_FAST : MAX_REF_CANDIDATES_PHASE2
                    if n_candidates > _max_k_greedy
                        overlap_scores = Vector{Tuple{Int,Int}}(undef, n_candidates)
                        for (ci2, rix) in enumerate(wstart:wend)
                            ov = _sorted_overlap_count(nl, all_nl[rix])
                            overlap_scores[ci2] = (ov, rix)
                        end
                        sort!(overlap_scores; by = x -> -x[1])
                        phase2_indices = [overlap_scores[k][2] for k in 1:min(_max_k_greedy, n_candidates)]
                    else
                        phase2_indices = collect(wstart:wend)
                    end

                    # Phase 2: full analytical evaluation on top candidates
                    for rix in phase2_indices
                        _merge_positions_adds!(nl, all_nl[rix], _pos_a, _adds_a)
                        bits, mil_val = _evaluate_candidate_greedy_analytical(_pos_a, _adds_a, params, T, _zz_vid, _mil_options)
                        if bits < best_bits
                            best_bits = bits
                            best_mil_val = mil_val
                            best_is_ref = true
                            best_ref_idx = rix
                        end
                    end
                end

                if best_is_ref
                    use_ref = true
                    ref_delta_val = UInt32(idx_local - best_ref_idx)
                    # Final merge for the winner — swap into storage
                    _merge_positions_adds!(nl, all_nl[best_ref_idx], _pos_a, _adds_a)
                    ref_positions_list[idx_local] = copy(_pos_a)
                    additions_list[idx_local] = copy(_adds_a)
                else
                    ref_positions_list[idx_local] = T[]
                    additions_list[idx_local] = T[]
                end
                mil_vec[idx_local] = best_mil_val
            elseif params.intra_ref_enabled && idx_local > 1
                # Analytical reference decision with 2-phase pruning + early termination
                _zz_vid = params.intra_zigzag ? T(idx_local) : nothing

                # Raw estimation (analytical)
                raw_bits, raw_use_iv = _estimate_raw_cost_analytical(nl, params, T, _zz_vid)
                best_raw_mode = raw_use_iv

                # Ref delta header overhead
                ref_overhead = 0
                if params.intra_ref_fixwidth
                    ref_overhead = max(1, ceil(Int, log2(params.intra_ref_window)))
                end

                wstart = max(1, idx_local - params.intra_ref_window)
                wend = idx_local - 1
                n_candidates = wend - wstart + 1

                # Phase 1: overlap screening
                _max_k_ref = params.cost_model == Compression.COST_MODEL_FAST ? MAX_REF_CANDIDATES_PHASE2_FAST : MAX_REF_CANDIDATES_PHASE2
                if n_candidates > _max_k_ref
                    overlap_scores = Vector{Tuple{Int,Int}}(undef, n_candidates)
                    for (ci2, rix) in enumerate(wstart:wend)
                        ov = _sorted_overlap_count(nl, all_nl[rix])
                        overlap_scores[ci2] = (ov, rix)
                    end
                    sort!(overlap_scores; by = x -> -x[1])
                    phase2_indices = [overlap_scores[k][2] for k in 1:min(_max_k_ref, n_candidates)]
                else
                    phase2_indices = collect(wstart:wend)
                end

                # Phase 2: analytical evaluation with early termination
                best_bits = raw_bits
                best_idx = 0
                # Use double-buffer swap: _pos_a/_adds_a for current, _pos_b/_adds_b for best
                for rix in phase2_indices
                    _merge_positions_adds!(nl, all_nl[rix], _pos_a, _adds_a)
                    ref_len = length(all_nl[rix])
                    bits, cm, aim = _evaluate_candidate_analytical(_pos_a, _adds_a, ref_len, params, T, _zz_vid; best_so_far=best_bits - ref_overhead)
                    total = bits + ref_overhead
                    if total < best_bits
                        best_bits = total
                        best_idx = rix
                        best_copy_mode = cm
                        best_add_mode = aim
                        # Swap buffers: _pos_b/_adds_b now hold the best result
                        _pos_a, _pos_b = _pos_b, _pos_a
                        _adds_a, _adds_b = _adds_b, _adds_a
                    end
                end
                if best_idx > 0
                    use_ref = true
                    ref_delta_val = UInt32(idx_local - best_idx)
                    # Best result is in _pos_b/_adds_b (after last swap)
                    ref_positions_list[idx_local] = copy(_pos_b)
                    additions_list[idx_local] = copy(_adds_b)
                end
                if !use_ref
                    ref_positions_list[idx_local] = T[]
                    additions_list[idx_local] = T[]
                end
            else
                ref_positions_list[idx_local] = T[]
                additions_list[idx_local] = T[]
            end
            use_ref_vec[idx_local] = use_ref
            ref_delta_vec[idx_local] = ref_delta_val
            copy_mode_vec[idx_local] = best_copy_mode
            add_mode_vec[idx_local] = best_add_mode
            raw_mode_vec[idx_local] = best_raw_mode
            if use_ref
                ref_index = idx_local - Int(ref_delta_val)
                ref_len_list[idx_local] = ref_index >= 1 ? length(all_nl[ref_index]) : 0
            else
                ref_len_list[idx_local] = 0
            end
        end


    end
    if parallel_search && s > 1 && Threads.nthreads(:default) > 1
        next = Threads.Atomic{Int}(1)
        @sync for _ in 1:min(search_workers, Threads.nthreads(:default), s)
            Threads.@spawn search_vertices!(next)
        end
        # Callbacks remain ordered and run on the caller, never worker tasks.
        if progress !== nothing
            for i in 1:s; progress(i, s, ci, n_clusters); end
        end
    else
        search_vertices!(nothing)
    end
    return (use_ref_vec, ref_delta_vec, ref_positions_list, additions_list, ref_len_list,
            copy_mode_vec, add_mode_vec, raw_mode_vec, mil_vec)
end
