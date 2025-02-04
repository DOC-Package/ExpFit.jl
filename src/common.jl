"""
    solve_vandermonde()


"""

function solve_vandermonde(
    hk::Vector{ComplexF64},
    gamm::Vector{ComplexF64},
    dt::Float64
)
    
    K = length(hk)
    M = length(gamm)
    
    # Solve the Vandermonde system
    G = Matrix{ComplexF64}(undef, K, M)
    G = [gamm[m]^k for k in 0:K-1, m in 1:M]

    # Coefficients
    coeff = G \ hk

    # Exponents
    exponent = Vector{ComplexF64}(undef, M)
    exponent = [-log(gamm[m]) / dt for m in 1:M]
    
    # Sort by magnitude
    idx = sortperm(coeff, by = x -> abs(x), rev=true)
    exponent = exponent[idx]
    coeff = coeff[idx]

    return exponent, coeff
end

"""
    hankel_matrix(hk; q=nothing)

Construct a Hankel matrix from the discrete data `hk` (a vector of complex numbers).

- `hk` : Data vector (of length K).
- `q`  : Number of columns for the Hankel matrix (optional). 
         If not provided, the default is:
           - For K = 2N (even), q is set to N+1.
           - For K = 2N-1 (odd),  q is set to N.
         The allowed values are:
           - For K = 2N: 2 ≤ q ≤ N+1.
           - For K = 2N-1: 2 ≤ q ≤ N.

Returns: A p × q Hankel matrix H such that H[i,j] = hk[i+j-1],
         where p is computed as p = K + 1 - q.
"""
function hankel_matrix(hk::AbstractVector{<:ComplexF64}; q::Union{Int,Nothing}=nothing) :: Matrix{ComplexF64}
    K = length(hk)
    
    if iseven(K)
        N = K ÷ 2
    else
        N = (K + 1) ÷ 2
    end

    # Set default q if not provided.
    if q === nothing
        if iseven(K)
            q = N + 1    
        else
            q = N    
        end
    end

    # Validate that q is within the allowed range.
    if iseven(K)
        @assert 2 <= q <= N+1 
    else
        @assert 2 <= q <= N 
    end

    # Compute the number of rows p.
    p = K + 1 - q

    return [hk[i + j - 1] for i in 1:p, j in 1:q]
end

function sumexp(t, a, c)
    s = 0.0 + 0.0im
    @inbounds for m in 1:length(a)
        s += c[m] * exp(-a[m] * t)
    end
    return s
end