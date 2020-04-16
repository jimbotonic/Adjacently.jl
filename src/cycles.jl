#
# JCNL: Julia Complex Networks Library
# Copyright (C) 2016-2020 Jimmy Dubuisson <jimmy.dubuisson@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#

using Graphs, DataStructures, Logging

"""
Get the list of colinks
"""
function list_colinks(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},rg::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},filename::String) where {T<:Unsigned}
	vs = vertices(g)
	n = length(vs)
	io = T[]
	# number of unique colinks per node
	un2 = zeros(Uint64,n)
	count = 0
	total = 0
	ttotal = 0
	bmax_size = 1000
	# delete file
	if isfile(filename)
		rm(filename)
	end
	for v in vs
		# get greater children and parents
		gc = filter(x->x>v,out_neighbors(v,g))
		gp = filter(x->x>v,out_neighbors(v,rg))
		inter = intersect(gc,gp)
		for t in inter
			append!(io,[v,t])
			ttotal += 1
		end
		un2[v] = length(inter)
		total += 1
		count += 1
		if count == bmax_size
			@debug("writing $bmax_size computed vertices [total # colinks: $ttotal | total # of computed vertices: $total (", (total/n*100),"%)]")
			f = open(filename, "a")
			write(f,io)
			io = T[]
			close(f)
			count = 0
		end
	end
	serialize_to_file(un2,"vertex-un2.jld")
	@debug("Total # of unique colinks: ", sum(un2))
	f = open(filename, "a")
	write(f,io)
	close(f)
end

"""
List triangles - map
"""
function list_triangles_map(v::T,g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},rg::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},dpos::Dict{T,T}) where {T<:Unsigned}
	candidates = (T,T,T)[]
	gc = filter(x->dpos[x]>dpos[v],out_neighbors(v,g))
	gp = filter(x->dpos[x]>dpos[v],out_neighbors(v,rg))
	#@debug("# greater children: ", length(gc))
	#@debug("# greater parents: ", length(gp))
	for c in gc
		# rank(child) > rank(v)
		for p in gp
			# child != parents && rank(parent) > rank(v)
			c != p && push!(candidates,(v,c,p))
		end
	end
	return candidates
end

"""
    list_triangles_reduce(candidates::Array{(T,T,T),1},g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},io::Array{T,1}) where {T<:Unsigned}

list triangles - reduce
"""
function list_triangles_reduce(candidates::Array{Tuple{T,T,T},1},g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},io::Array{T,1}) where {T<:Unsigned}
	count = convert(T,0)
	for c in candidates
		cc = out_neighbors(convert(T,c[2]),g)
		sort!(cc)
		if binary_search(cc,c[3]) != -1
			append!(io,[c[1],c[2],c[3]])
			count = convert(T,1+count)
		end
	end
	return count
end

"""
List triangles - map&reduce
"""
function list_triangles_mapreduce(v::T,g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},rg::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},dpos::Dict{T,T},io::Array{T,1}) where {T<:Unsigned}
	count = convert(T,0)
	vpos = dpos[v]
	gc = filter(x->dpos[x]>vpos,out_neighbors(v,g))
	gp = filter(x->dpos[x]>vpos,out_neighbors(v,rg))
	for c in gc
		cc = out_neighbors(c,g)
		sort!(cc)
		for p in gp
			if c != p
				# NB: searchsortedfirst is unreliable!
				#if searchsortedfirst(cc,p) <= length(cc)
				if binary_search(cc,p) != -1
					append!(io,[v,c,p])
					count = convert(T,1+count)
				end
			end
		end
	end
	return count
end

"""
Initialize un2 (# of colink candidates) and un3 (# of triangle candidates) array
"""
function init_un23_arrays(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},rg::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},dpos::Dict{T,T}) where {T<:Unsigned}
	@debug("initializing # colink candidates + # triangle candidates arrays")
	vs = vertices(g)
	n = length(vs)
	un2 = Uint64[]
	un3 = Uint64[]
	for v in vs
		# total order relation 1
		# get the set of greater children of vertex v
		gc1 = filter(x->x>v,out_neighbors(v,g))
		push!(un2, length(gc1))

		# total order relation 2
		vpos = dpos[v]
		gc2 = filter(x->dpos[x]>vpos,out_neighbors(v,g))
		gp = filter(x->dpos[x]>vpos,out_neighbors(v,rg))
		p = length(gp)
		c = length(gc2)
		i = length(intersect(gp,gc2))
		push!(un3, (p*c-i))
	end
	return un2,un3
end

"""
Initialize mn2 (maximum # of colinks) and mn3 (maximum # of triangles) array
"""
function init_mn23_arrays(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},rg::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}}) where {T<:Unsigned}
	@debug("initializing max # colinks + max # triangles arrays")
	vs = vertices(g)
	n = length(vs)
	mn2 = Uint64[]
	mn3 = Uint64[]
	for v in vs
		parents = out_neighbors(v,rg)
		children = out_neighbors(v,g)
		p = length(parents)
		c = length(children)
		i = length(intersect(parents,children))
		# maximum number of triangles:
		# = (p_minus_c * c_minus_p) + i * (p_minus_c + c_minus_p + (i - 1))
		# = p*c - i
		push!(mn2, (p+c-i))
		push!(mn3, (p*c-i))
	end
	return mn2,mn3
end

"""
Get the list of triangles
"""
function list_triangles(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},rg::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},filename::String) where {T<:Unsigned}
	n = length(vertices(g))
	io = T[]
	# actual number of unique triangles
	un3 = zeros(Uint64,n)
	count = 0
	total = 0
	ttotal = 0
	bmax_size = 10000000
	# delete triangle file if needed
	if isfile(filename)
		rm(filename)
	end
	# maximum number of triangles
	#mn3 = load_serialized("vertex-mn3.jld")
	mn2,mn3 = init_mn23_arrays(g,rg)
	#serialize_to_file(mn2,"vertex-mn2.jld")
	#serialize_to_file(mn3,"vertex-mn3.jld")
	@debug("sorting mn3 array")
	#smn3 = load_serialized("smn3.jld")
	smn3 = bottom_up_sort!(mn3)
	#serialize_to_file(smn3,"smn3.jld")
	@debug("mn3 array sorted")
	@debug("loading position dictionary")
	#dpos = load_serialized("dpos.jld")
	dpos = Dict{T,T}()
	pos = 1
	for v in smn3
		dpos[v] = pos
		pos += 1
	end
	#serialize_to_file(dpos,"dpos.jld")
	@debug("dictionary loaded")
	for v in smn3
		#@debug("Checking vertex: $v")
		#tic()
		#candidates = list_triangles_map(convert(T,v),g,rg,dpos)
		#@debug("# candidates: ", length(candidates))
		#nadded = list_triangles_reduce(candidates,g,io)
		nadded = list_triangles_mapreduce(convert(T,v),g,rg,dpos,io)
		un3[v] = nadded
		#@debug("# found triangles for vertex $v: $nadded")
		#@debug(toc())
		ttotal += nadded
		count += nadded
		total += 1
		if count >= bmax_size
			@debug("writing $count triangles [total # triangles: $ttotal | total # of computed vertices: $total (", (total/n*100),"%)]")
			f = open(filename, "a")
			write(f,io)
			io = T[]
			close(f)
			count = 0
		end
	end
	serialize_to_file(un3,"vertex-un3.jld")
	@debug("Total # of unique triangles: ", sum(un3))
	# write remaining triangles
	f = open(filename, "a")
	write(f,io)
	close(f)
end

# load the colinks distribution from the specified file
function load_colinks_distribution(size::T,filename::String) where {T<:Unsigned}
	f = open(filename,"r")
	d = Dict{T,T}()
	while !eof(f)
		v1 = reinterpret(T,read(f,Uint8,sizeof(T)))[1]
		v2 = reinterpret(T,read(f,Uint8,sizeof(T)))[1]
		if haskey(d,v1)
			d[v1] += 1
		else
			d[v1] = 1
		end
		if haskey(d,v2)
			d[v2] += 1
		else
			d[v2] = 1
		end
	end
	close(f)
	c_stats = zeros(T,size)
	for i in 1:size
		if haskey(d,i)
			c_stats[i] = d[i]
		else
			c_stats[i] = 0
		end
	end
	return c_stats
end

# load the triangles distribution from the specified file
function load_triangles_distribution(size::T,filename::String) where {T<:Unsigned}
	f = open(filename,"r")
	d = Dict{T,T}()
	while !eof(f)
		v1 = reinterpret(T,read(f,Uint8,sizeof(T)))[1]
		v2 = reinterpret(T,read(f,Uint8,sizeof(T)))[1]
		v3 = reinterpret(T,read(f,Uint8,sizeof(T)))[1]
		if haskey(d,v1)
			d[v1] += 1
		else
			d[v1] = 1
		end
		if haskey(d,v2)
			d[v2] += 1
		else
			d[v2] = 1
		end
		if haskey(d,v3)
			d[v3] += 1
		else
			d[v3] = 1
		end
	end
	close(f)
	t_stats = zeros(T,size)
	for i in 1:size
		if haskey(d,i)
			t_stats[i] = d[i]
		else
			t_stats[i] = 0
		end
	end
	return t_stats
end
