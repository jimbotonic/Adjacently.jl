#
# Zeta encoding on Amazon graph in children and index modes,
# with and without reference encoding; verify decoded children match.
#

include("run_tests_main.jl")

@testset "Amazon zeta: children/index with and without reference" begin
    @info "=== Loading Amazon dataset for zeta children/index tests ==="

    # Load full Amazon graph, then take a manageable, properly remapped subset
    amz_g = load_dataset(AMZ_DATASET_IN; separator='\t')
    amz_core, _, _ = get_core(amz_g)
    full_neighbor_lists = get_neighbor_lists(amz_core)

    # Use compact integer type and zeta encoding
    T = UInt24
    encoding = :zeta

    # Create a small, contiguous subset to keep runtime reasonable
    vertex_ids = sort(collect(keys(full_neighbor_lists)))
    subset_size = length(vertex_ids)
    selected_vertex_ids = vertex_ids

    # Remap vertex IDs 1..subset_size and filter neighbors to the subset
    old_to_new = Dict{T,T}()
    for (new_id, old_id) in enumerate(selected_vertex_ids)
        old_to_new[old_id] = T(new_id)
    end

    remapped_neighbor_lists = Dict{T, Vector{T}}()
    for (new_id, old_id) in enumerate(selected_vertex_ids)
        old_neighbors = full_neighbor_lists[old_id]
        new_neighbors = T[]
        for old_neighbor in old_neighbors
            if haskey(old_to_new, old_neighbor)
                push!(new_neighbors, old_to_new[old_neighbor])
            end
        end
        sort!(new_neighbors)
        remapped_neighbor_lists[T(new_id)] = new_neighbors
    end

    total_edges = sum(length(neigh) for neigh in values(remapped_neighbor_lists))
    @info "Subset: vertices=$(subset_size), edges=$(total_edges)"

    for mode in (:children, :index)
        for use_reference in (false, true)
            @testset "zeta-$mode-ref=$(use_reference)" begin
                @info "Encoding: mode=$mode, reference=$(use_reference)"

                # Encode
                io = IOBuffer()
                writer = BitWriter(io)
                Adjacently.Compression.write_compressed_graph_data(writer, remapped_neighbor_lists, encoding, mode, use_reference)
                flush_bitwriter(writer; flush_last_bits=true)

                bytes = position(io)
                bpe = total_edges > 0 ? (bytes * 8) / total_edges : 0.0
                @info "  Size: $(bytes) bytes ($(round(bpe, digits=3)) bits/edge)"

                # Decode
                seekstart(io)
                reader = BitReader(io)
                decoded = Adjacently.Compression.read_compressed_graph_data(reader, T(subset_size), encoding, mode, T)

                # Verify: compare children for every vertex
                mismatches = 0
                for v in 1:subset_size
                    vid = T(v)
                    expected = remapped_neighbor_lists[vid]
                    actual = get(decoded, vid, T[])
                    if Set(expected) != Set(actual)
                        mismatches += 1
                        if mismatches <= 3
                            @error "  Mismatch v=$v: expected $(expected), got $(actual)"
                        end
                    end
                end

                if mismatches == 0
                    @info "  ✅ All $(subset_size) vertices match"
                    @test true
                else
                    @error "  ❌ $mismatches / $(subset_size) vertices mismatched"
                    @test false
                end
            end
        end
    end
end
