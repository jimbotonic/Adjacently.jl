#
# Graph Attention Network (GAT) layer with ultimate numerical stability.
#

using LinearAlgebra, SparseArrays, Random

mutable struct GATLayer
    W::Matrix{Float64}      # (in_dim, out_dim)
    a::Matrix{Float64}      # (2 * out_dim, 1)
    in_dim::Int
    out_dim::Int
    alpha::Float64          # LeakyReLU slope
    dropout::Float64        # Dropout rate
    concat::Bool            # Concatenate heads vs. averaging
    # Caches
    Z::Matrix{Float64}      # (N, out_dim)
    e_unactivated::Matrix{Float64} # (num_edges, 1)
    att::SparseMatrixCSC{Float64, Int} # Attention weights
    # Gradients
    dW::Matrix{Float64}
    da::Matrix{Float64}
end

function GATLayer(in_dim::Int, out_dim::Int; alpha=0.2, dropout=0.5, concat=true, seed=42)
    rng = MersenneTwister(seed)
    W = randn(rng, in_dim, out_dim) .* sqrt(2.0 / (in_dim + out_dim))
    a = randn(rng, 2 * out_dim, 1) .* sqrt(2.0 / (2 * out_dim + 1))
    return GATLayer(W, a, in_dim, out_dim, alpha, dropout, concat,
                    zeros(0,0), zeros(0,0), spzeros(0,0), zeros(size(W)), zeros(size(a)))
end

function leaky_relu(x, alpha)
    return max.(alpha .* x, x)
end

function leaky_relu_derivative(x, alpha)
    return Float64.(x .> 0) .+ alpha .* Float64.(x .<= 0)
end

function gat_forward(layer::GATLayer, X::Matrix{Float64}, adj::SparseMatrixCSC; training::Bool=true)
    N = size(X, 1)
    layer.Z = X * layer.W
    if any(isnan, layer.Z); layer.Z = replace(layer.Z, NaN => 0.0); end

    rows, cols, _ = findnz(adj)
    e_unactivated = zeros(length(rows), 1)
    for idx in 1:length(rows)
        logit = dot(view(layer.a, 1:layer.out_dim), view(layer.Z, rows[idx], :)) +
                dot(view(layer.a, layer.out_dim+1:2*layer.out_dim), view(layer.Z, cols[idx], :))
        e_unactivated[idx] = logit
    end
    
    clamp!(e_unactivated, -10.0, 10.0)
    layer.e_unactivated = e_unactivated
    
    e_vals = exp.(leaky_relu(e_unactivated, layer.alpha))
    exp_adj = sparse(rows, cols, vec(e_vals), N, N)
    row_sums = exp_adj * ones(N)
    
    layer.att = copy(exp_adj)
    for j in 1:N
        for idx in layer.att.colptr[j]:(layer.att.colptr[j+1]-1)
            i = layer.att.rowval[idx]
            if row_sums[i] > 1e-12
                layer.att.nzval[idx] /= row_sums[i]
            else
                layer.att.nzval[idx] = 1e-6
            end
        end
    end

    if training && layer.dropout > 0.0
        mask = rand(size(layer.att.nzval)) .> layer.dropout
        layer.att.nzval .*= mask
    end

    H_prime = layer.att * layer.Z
    return layer.concat ? H_prime : leaky_relu(H_prime, layer.alpha)
end

function gat_backward!(layer::GATLayer, dH_prime, X, adj)
    N = size(X, 1)
    # Aggressive NaN protection
    if any(isnan, dH_prime); dH_prime = replace(dH_prime, NaN => 0.0); end
    clamp!(dH_prime, -5.0, 5.0)
    
    if !layer.concat
        dH_prime = dH_prime .* leaky_relu_derivative(layer.att * layer.Z, layer.alpha)
    end
    
    dZ = layer.att' * dH_prime
    
    row_inner_prods = zeros(N)
    for j in 1:N
        Zj = view(layer.Z, j, :)
        for idx in layer.att.colptr[j]:(layer.att.colptr[j+1]-1)
            i = layer.att.rowval[idx]
            datt_ij = dot(view(dH_prime, i, :), Zj)
            row_inner_prods[i] += datt_ij * layer.att.nzval[idx]
        end
    end
    
    datt_unsoftmax_vals = zeros(length(layer.att.nzval))
    for j in 1:N
        Zj = view(layer.Z, j, :)
        for idx in layer.att.colptr[j]:(layer.att.colptr[j+1]-1)
            i = layer.att.rowval[idx]
            datt_ij = dot(view(dH_prime, i, :), Zj)
            datt_unsoftmax_vals[idx] = layer.att.nzval[idx] * (datt_ij - row_inner_prods[i])
        end
    end
    
    # Clip attention gradients
    clamp!(datt_unsoftmax_vals, -5.0, 5.0)
    
    de_vals = datt_unsoftmax_vals .* leaky_relu_derivative.(layer.e_unactivated, layer.alpha)
    for i in 1:length(de_vals)
        if layer.e_unactivated[i] < -10.0 || layer.e_unactivated[i] > 10.0
            de_vals[i] = 0.0
        end
    end

    rows, cols, _ = findnz(adj)
    layer.da .= 0.0
    dZ_att = zeros(size(layer.Z))
    for idx in 1:length(rows)
        de = de_vals[idx]
        ri, ci = rows[idx], cols[idx]
        for d in 1:layer.out_dim
            layer.da[d] += de * layer.Z[ri, d]
            layer.da[layer.out_dim + d] += de * layer.Z[ci, d]
            dZ_att[ri, d] += de * layer.a[d]
            dZ_att[ci, d] += de * layer.a[layer.out_dim + d]
        end
    end
    dZ .+= dZ_att

    layer.dW = X' * dZ
    if any(isnan, layer.dW); layer.dW = replace(layer.dW, NaN => 0.0); end
    if any(isnan, layer.da); layer.da = replace(layer.da, NaN => 0.0); end
    max_grad = 5.0
    nw = norm(layer.dW); if nw > max_grad; layer.dW .*= max_grad / nw; end
    na = norm(layer.da); if na > max_grad; layer.da .*= max_grad / na; end

    return dZ * layer.W'
end
