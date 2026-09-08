# Usage: julia --project=. --threads=1 bench/graph_compression/nonra_decode.jl cnr_cg path/to/file.mgz
using Adjacently, JSON, SHA, Random, Statistics, Logging
global_logger(ConsoleLogger(stderr, Logging.Warn))
root=normpath(joinpath(@__DIR__,"../.."))
cfg=JSON.parsefile(joinpath(@__DIR__,"configs/nonra_adjacency_v1.json"))
length(ARGS)==2 || error("usage: nonra_decode.jl <eat_bg|eat_cs|eat_cg|cnr_bg|cnr_cs|cnr_cg> <file.mgz>")
key,path=ARGS
if key=="cnr_cg_best"
    best=JSON.parsefile(joinpath(@__DIR__,"configs/cg_cnr_nonra_v2.json"))
    cfg["files_sha256"][key]=best["sha256"]
end
Threads.nthreads(:default)==cfg["decode_workers"] || error("run with --threads=1")
haskey(cfg["files_sha256"],key) || error("unknown pinned file: $key")
digest(p)=bytes2hex(sha256(read(p)))
digest(path)==cfg["files_sha256"][key] || error("file does not match the pinned profile")
for line in eachline(joinpath(@__DIR__,"configs/nonra_adjacency_v1.sha256"))
    hash,relative=split(line;limit=2)
    digest(joinpath(root,strip(relative)))==hash || error("source drift: $(strip(relative)); use a clean matching checkout")
end
M=Adjacently.MGS
options=cfg["recommended"][occursin("_cg",key) ? "cg" : "bg_cs"]
flags=(;[Symbol(k)=>v for (k,v) in options if k!="compact_output"]...)
baseline=M.load_indexed_mgs3_graph(path;flags...)
variants=Dict(
    "legacy_preload"=>()->M.load_indexed_mgs3_graph(path;flags...),
    "adjacency"=>()->M.load_adjacency_mgs3_graph(path;flags...),
    "compact_output"=>()->M.load_adjacency_mgs3_graph(path;flags...,compact_output=true))
function verify(value)
    value isa AbstractVector || return
    @assert length(value)==baseline.n
    @assert all(value[v]==baseline.neighbors(v) for v in 1:baseline.n)
end
for fn in values(variants);verify(fn());end
samples=Any[];rng=MersenneTwister(53)
for rep in 1:7,variant in shuffle(rng,sort(collect(keys(variants))))
    GC.gc();t=@timed variants[variant]();verify(t.value)
    push!(samples,(;variant,rep,seconds=t.time,allocated_bytes=t.bytes,gc_seconds=t.gctime))
end
retained=Dict(k=>Base.summarysize(variants[k]()) for k in ("adjacency","compact_output"))
@assert digest(path)==cfg["files_sha256"][key]
JSON.print(stdout,(;profile=cfg["id"],file=key,sha256=digest(path),julia=string(VERSION),
    threads=Threads.nthreads(:default),n=baseline.n,m=baseline.m,bytes=filesize(path),
    bpe=8filesize(path)/baseline.m,flags,samples,retained_bytes_estimate=retained,
    query_semantics="Full-file preprocessing, not compressed random access"),2)
println()
