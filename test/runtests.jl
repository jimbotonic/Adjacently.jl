using Test

module CodecReaderTests
include("run_tests_codec_readers.jl")
end

# Independent modules prevent fixture/helper-name collisions between suites.
# JULIA_NUM_THREADS=4 julia --project -e 'using Pkg; Pkg.test()'
@testset "Opt-in codec acceleration" begin
    @testset "BG/CS exact streams" begin
        include("run_tests_parallel_encoding.jl")
    end
end
module CGOptionsTests
include("run_tests_cg_decode_options.jl")
include("run_tests_codec_promotion.jl")
end
module LLPOptionsTests
include("run_tests_parallel_llp.jl")
end
module AdjacencyLoaderTests
include("run_tests_adjacency_loader.jl")
end
