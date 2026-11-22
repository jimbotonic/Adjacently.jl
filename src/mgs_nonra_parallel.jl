# Opt-in non-RA two-pass reconstruction. Entropy parsing is serial; each
# task owns a forest traversal and scratch. Staging adds O(n+m) live memory.
# No original graph, decoded cache, checkpoint, or on-disk change is involved.

const _NONRA_ANCHOR=UInt64(1)<<40
struct _NonRACSParsed
    parents::Vector{Int}
    degrees::Vector{Int}
    copy_offsets::Vector{Int}
    copies::Vector{UInt32}
    residual_offsets::Vector{Int}
    coordinates::Vector{UInt64}
    modes::Dict{String,Int}
end
function _nonra_cs_hybrid!(a,count,rc,flag)
    empty!(a);hybrid=Compression.brc_decode_bit!(flag)
    first=UInt64(Int64(_NONRA_ANCHOR)+Compression._zigzag_decode(Compression.rc_decode_value!(rc)-1));push!(a,first)
    if !hybrid
        for _ in 2:count;push!(a,a[end]+Compression.rc_decode_value!(rc));end
    else
        sections=Int(Compression.rc_decode_value!(rc))
        for _ in 1:sections
            if !Compression.brc_decode_bit!(flag)
                n=Int(Compression.rc_decode_value!(rc))
                for _ in 1:n
                    length(a)>=count && break
                    push!(a,a[end]+Compression.rc_decode_value!(rc))
                end
            elseif !Compression.brc_decode_bit!(flag)
                n=Int(Compression.rc_decode_value!(rc))
                for _ in 1:n
                    value=Compression.rc_decode_value!(rc);len=Int(Compression.rc_decode_value!(rc))
                    for _ in 1:len
                        length(a)>=count && break
                        push!(a,a[end]+value)
                    end
                end
            else
                n=Int(Compression.rc_decode_value!(rc))
                for _ in 1:n
                    start=Compression.rc_decode_value!(rc);len=Int(Compression.rc_decode_value!(rc))
                    length(a)>=count && break
                    # Hybrid interval starts are absolute rank-space coordinates,
                    # not offsets from the vertex. Keep them unshifted.
                    push!(a,start)
                    for j in 1:len-1
                        length(a)>=count && break
                        push!(a,start+UInt64(j))
                    end
                end
            end
            length(a)>=count && break
        end
    end
    @assert length(a)==count
    a
end
function _nonra_cs_parse(f)
    (f.alg == :cs && f.p.lr_split) || throw(ArgumentError("CS staging requires LR split"))
    0 < f.n <= typemax(UInt32) || throw(ArgumentError("CS staging supports at most 2^32-1 vertices"))
    r=BitReader(f.structural)
    Compression.read_bit(r) && throw(ArgumentError("CS staging does not support indexed streams"))
    ie=get(Compression.TAG_ENCODINGS,Compression._read_encoding_tag(r),:fibonacci)
    rc=Compression.CtxRangeDecoder(f.resid);rd=Compression.CtxRangeDecoder(f.rd)
    cp=Compression.BinRangeDecoder(f.cp);cmd=Compression.CtxRangeDecoder(f.cmd);flag=Compression.BinRangeDecoder(f.flag)
    n=f.n;parents=zeros(Int,n);degrees=zeros(Int,n)
    co=Vector{Int}(undef,n+1);ro=Vector{Int}(undef,n+1)
    copies=UInt32[];coords=UInt64[];scratch=UInt64[];merge=UInt64[];modecounts=zeros(Int,7)
    for v in 1:n
        co[v]=length(copies)+1;ro[v]=length(coords)+1
        is_empty,refmode,encoding,mil=Compression._cs_cmd_decode(Int(Compression.rc_decode_value!(cmd)))
        if is_empty;modecounts[1]+=1;continue;end
        ei=encoding==:interval ? 1 : encoding==:delta ? 2 : encoding==:rle ? 3 : error("Unsupported grammar")
        modecounts[1+ei+(refmode==:none ? 0 : 3)]+=1
        if refmode!=:none
            p=v-Int(Compression.rc_decode_value!(rd));@assert max(1,v-f.p.ref_window_size)<=p<v
            parents[v]=p
            for pos in 1:degrees[p]
                Compression.brc_decode_bit!(cp) && push!(copies,UInt32(pos))
            end
        end
        Compression.rc_reset_region!(rc)
        residuals=if encoding==:interval
            Compression._read_intervals_lr_source(r,ie,mil,UInt64,_NONRA_ANCHOR,f.p.tight_intervals,rc,scratch,merge)
        elseif encoding==:delta
            Compression._read_nonra_stop(r,ie,UInt64,_NONRA_ANCHOR,rc,flag,scratch)
        elseif encoding==:rle
            count=Int(Compression.rc_decode_value!(rc))-1
            count==0 ? empty!(scratch) : _nonra_cs_hybrid!(scratch,count,rc,flag)
        else
            error("Unsupported residual grammar $encoding")
        end
        append!(coords,residuals)
        degrees[v]=length(copies)+1-co[v]+length(residuals)
    end
    co[n+1]=length(copies)+1;ro[n+1]=length(coords)+1
    labels=("empty","none/interval","none/delta","none/rle","reference/interval","reference/delta","reference/rle")
    modes=Dict(labels[i]=>modecounts[i] for i in 1:7)
    _NonRACSParsed(parents,degrees,co,copies,ro,coords,modes)
end
function _nonra_forest(t)
    n=length(t.parents);counts=zeros(Int,n);roots=Int[]
    for (v,p) in enumerate(t.parents);p==0 ? push!(roots,v) : (counts[p]+=1);end
    offsets=Vector{Int}(undef,n+1);offsets[1]=1
    for v in 1:n;offsets[v+1]=offsets[v]+counts[v];end
    cursor=copy(offsets);children=Vector{Int}(undef,offsets[end]-1)
    for (v,p) in enumerate(t.parents)
        if p>0;children[cursor[p]]=v;cursor[p]+=1;end
    end
    (;roots,offsets,children)
end
function _nonra_cs_vertex!(adj::Vector{Vector{T}},t,v,scratch) where T
    p=t.parents[v];copied=T[]
    if p>0
        ref=adj[p];sizehint!(copied,t.copy_offsets[v+1]-t.copy_offsets[v])
        for i in t.copy_offsets[v]:t.copy_offsets[v+1]-1;push!(copied,ref[Int(t.copies[i])]);end
    end
    vid=isempty(copied) ? T(v) : Compression._rank_of(T(v),copied)
    adds=p==0 ? T[] : empty!(scratch)
    for i in t.residual_offsets[v]:t.residual_offsets[v+1]-1
        x=t.coordinates[i]
        value=x>_NONRA_ANCHOR÷2 ? Int64(x)-Int64(_NONRA_ANCHOR)+Int64(vid) : Int64(x)
        push!(adds,T(value))
    end
    if p>0
        !isempty(copied) && Compression._rankgap_inverse!(adds,copied)
        adj[v]=Compression._merge_sorted_owned!(copied,adds)
    else
        adj[v]=adds
    end
    @assert length(adj[v])==t.degrees[v]
end
function _nonra_cs_reconstruct(t,f,workers,::Type{T}) where T
    adj=Vector{Vector{T}}(undef,length(t.parents));next=Threads.Atomic{Int}(1)
    function worker()
        stack=Int[];scratch=T[]
        while true
            start=Threads.atomic_add!(next,32);start>length(f.roots) && break
            for i in start:min(start+31,length(f.roots))
                push!(stack,f.roots[i])
                while !isempty(stack)
                    v=pop!(stack);_nonra_cs_vertex!(adj,t,v,scratch)
                    for j in f.offsets[v]:f.offsets[v+1]-1;push!(stack,f.children[j]);end
                end
            end
        end
    end
    if workers==1;worker()
    else;@sync for _ in 1:min(workers,Threads.nthreads(:default));Threads.@spawn worker();end
    end
    adj
end

struct _NonRACGTokens
    parents::Vector{Int}
    degrees::Vector{Int}
    copy_offsets::Vector{Int}
    copies::Vector{UInt32}
    token_offsets::Vector{Int}
    tokens::Vector{UInt64}
    mil::Int
end
@inline function _nonra_token!(values,decoder)
    v=Compression.rc_decode_value!(decoder);push!(values,v);v
end
function _nonra_cg_parse(f)
    p=f.p;n=f.n
    (0 < n <= typemax(UInt32) && p.membership == :implicit_ranges &&
        p.intra_lr_split && p.intra_intervals && p.intra_zigzag &&
        !p.intra_tight_deltas && p.intra_ref_enabled && p.intra_ref_fixwidth &&
        !p.intra_ref_vlc && !p.intra_greedy_mil && !p.intra_block_try && !p.intra_mgs) ||
        throw(ArgumentError("Unsupported CG staged reconstruction parameters"))
    r=BitReader(f.structural)
    clusters=Compression.CG._read_membership(r,p,UInt32)
    (length(clusters)==1 && length(clusters[1])==n &&
        all(i -> clusters[1][i] == i, 1:n)) ||
        throw(ArgumentError("CG staging requires identity K=1 membership"))
    rd=f.rd;cp=f.cp;resid=f.resid
        refdecoder=Compression.CtxRangeDecoder(rd);copydecoder=Compression.BinRangeDecoder(cp);decoder=Compression.CtxRangeDecoder(resid)
        used,deltas=Compression.CG._read_fixwidth_ref_deltas(r,n,p.intra_ref_window;vlc=false,refdist_source=refdecoder)
        parents=zeros(Int,n);degrees=zeros(Int,n)
        co=Vector{Int}(undef,n+1);to=Vector{Int}(undef,n+1)
        copies=UInt32[];tokens=UInt64[]
        for v in 1:n
            co[v]=length(copies)+1;to[v]=length(tokens)+1
            if used[v]
                parent=v-Int(deltas[v]);@assert 1<=parent<v;parents[v]=parent
                for pos in 1:degrees[parent]
                    Compression.brc_decode_bit!(copydecoder) && push!(copies,UInt32(pos))
                end
            end
            degree=length(copies)+1-co[v]
            Compression.rc_reset_region!(decoder)
            ni=Int(_nonra_token!(tokens,decoder))-1
            for _ in 1:ni
                _nonra_token!(tokens,decoder)
                degree+=Int(_nonra_token!(tokens,decoder))-1+p.intra_mil
            end
            nr=Int(_nonra_token!(tokens,decoder))-1;degree+=nr
            if nr>0
                left=Int(_nonra_token!(tokens,decoder))-1;@assert 0<=left<=nr
                for _ in 1:nr;_nonra_token!(tokens,decoder);end
            end
            degrees[v]=degree
        end
        co[n+1]=length(copies)+1;to[n+1]=length(tokens)+1
        @assert isempty(Compression.CG.read_stop_delta_list(r;encoding=p.gap,T=UInt32))
        _NonRACGTokens(parents,degrees,co,copies,to,tokens,p.intra_mil)
end
function _nonra_cg_vertex!(adj::Vector{Vector{T}},t::_NonRACGTokens,v,add_scratch,merge_scratch) where T
    p=t.parents[v];vid=T(v)
    copied=T[]
    if p>0
        ref=adj[p]
        sizehint!(copied,t.copy_offsets[v+1]-t.copy_offsets[v])
        for i in t.copy_offsets[v]:t.copy_offsets[v+1]-1
            push!(copied,ref[Int(t.copies[i])])
        end
        !isempty(copied) && (vid=T(Compression._rank_of(vid,copied)))
    end
    adds=p==0 ? T[] : empty!(add_scratch)
    values=t.tokens;i=t.token_offsets[v]
    ni=Int(values[i])-1;i+=1;previous=zero(T)
    for j in 1:ni
        start=j==1 ? T(Int64(vid)+Compression._zigzag_decode(values[i]-1)) : previous+T(values[i]);i+=1
        len=Int(values[i])-1+t.mil;i+=1
        for k in 0:len-1;push!(adds,start+T(k));end
        previous=start
    end
    split=length(adds);nr=Int(values[i])-1;i+=1
    if nr>0
        left=Int(values[i])-1;i+=1
        sizehint!(adds,split+nr;shrink=false)
        distance=zero(T)
        for j in 1:left
            distance=j==1 ? T(values[i]) : distance+T(values[i])-one(T);i+=1
            push!(adds,vid-distance)
        end
        reverse!(adds,split+1,length(adds))
        distance=zero(T)
        for j in 1:nr-left
            distance=j==1 ? T(values[i]) : distance+T(values[i])-one(T);i+=1
            push!(adds,vid+distance-one(T))
        end
    end
    @assert i==t.token_offsets[v+1]
    Compression._merge_lr_runs!(adds,split,merge_scratch)
    if p>0
        !isempty(copied) && Compression._rankgap_inverse!(adds,copied)
        adj[v]=Compression._merge_sorted_owned!(copied,adds)
    else
        adj[v]=adds
    end
    @assert length(adj[v])==t.degrees[v]
    nothing
end

function _nonra_cg_reconstruct(t::_NonRACGTokens,f,workers,::Type{T}) where T
    n=length(t.parents);adj=Vector{Vector{T}}(undef,n)
    next=Threads.Atomic{Int}(1)
    function worker()
        stack=Int[];adds=T[];merge=T[]
        while true
            firstroot=Threads.atomic_add!(next,32)
            firstroot>length(f.roots) && break
            for rootindex in firstroot:min(firstroot+31,length(f.roots))
                push!(stack,f.roots[rootindex])
                while !isempty(stack)
                    v=pop!(stack);_nonra_cg_vertex!(adj,t,v,adds,merge)
                    for i in f.offsets[v]:f.offsets[v+1]-1;push!(stack,f.children[i]);end
                end
            end
        end
    end
    if workers==1
        worker()
    else
        @sync for _ in 1:min(workers,Threads.nthreads(:default));Threads.@spawn worker();end
    end
    adj
end


function _nonra_parallel_adjacency(f, workers::Int)
    workers > 0 || throw(ArgumentError("decode_workers must be positive"))
    if f.alg == :cs
        parsed = _nonra_cs_parse(f)
        tree = _nonra_forest(parsed)
        return _nonra_cs_reconstruct(parsed, tree, workers, f.T)
    elseif f.alg == :cg
        parsed = _nonra_cg_parse(f)
        tree = _nonra_forest(parsed)
        return _nonra_cg_reconstruct(parsed, tree, workers, f.T)
    end
    throw(ArgumentError("Parallel reconstruction supports CS and CG only"))
end
