module Constants

export BUFFER_SIZE, FIB_NUMBERS, GOLOMB_BASE, ZETA_BASE, ZETA_H_BOUNDS, ZETA_POWER_BASES

const BUFFER_SIZE::Int = 4096
const GOLOMB_BASE::Int = 64
const ZETA_BASE::Int = 4

# Precomputed Fibonacci numbers <= 2^40
# NB: this is used for the Fibonacci code
const FIB_NUMBERS::Vector{UInt} = let
    fib = UInt[1, 2]
    while true
        next = fib[end] + fib[end - 1]
        next > UInt(2)^40 && break
        push!(fib, next)
    end
    fib
end

# Precomputed zeta bounds for fast h calculation
# ZETA_H_BOUNDS[k] = vector of threshold values where h increments
# For h=0: v in [1, 2^k-1], for h=1: v in [2^k, 2^(2k)-1], etc.
# ZETA_H_BOUNDS[k][h+1] = 2^(k*(h+1)) = first value that requires h+1
const ZETA_H_BOUNDS = let
    max_k = 8  # Support more k values since we're only storing bounds
    max_bits = 40  # Support up to 2^40 as requested
    
    bounds = Dict{Int, Vector{UInt64}}()
    
    for k in 1:max_k
        bound_list = UInt64[]
        h = 0
        
        while true
            k_times_h_plus_1 = k * (h + 1)
            if k_times_h_plus_1 > max_bits
                break
            end
            
            # 2^(k*(h+1)) = threshold where h transitions to h+1
            threshold = UInt64(1) << k_times_h_plus_1
            push!(bound_list, threshold)
            h += 1
            
            # Stop if we've reached the maximum value we can represent
            if k_times_h_plus_1 >= 63  # UInt64 limit
                break
            end
        end
        
        bounds[k] = bound_list
    end
    
    bounds
end

# Precomputed power_base values for zeta coding  
# ZETA_POWER_BASES[k][h+1] = 2^(k*h) for k and h (1-based indexing)
const ZETA_POWER_BASES = let
    max_k = 8  # Match ZETA_H_BOUNDS
    max_h = 64  # Maximum h value we might encounter
    
    bases = Dict{Int, Vector{UInt64}}()
    
    for k in 1:max_k
        power_list = UInt64[]
        
        for h in 0:max_h
            k_times_h = k * h
            if k_times_h <= 63  # Avoid UInt64 overflow
                push!(power_list, UInt64(1) << k_times_h)
            else
                break
            end
        end
        
        bases[k] = power_list
    end
    
    bases
end

end # module Constants 