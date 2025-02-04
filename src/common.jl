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
    hankel_matrix(hk; p=nothing)

Construct a Hankel matrix from the discrete data `hk` (a vector of complex numbers).

- `hk` : Data vector (of length K).
- `p`  : Number of rows for the Hankel matrix (optional).
         If not provided, `p` is set to ⌊(K+1)/2⌋ and `q` is computed as `K+1-p`.

Returns: A p × q Hankel matrix H such that H[i,j] = hk[i+j-1].
"""
function hankel_matrix(hk::AbstractVector{<:ComplexF64}; p::Union{Int,Nothing}=nothing)
    K = length(hk)
    if p === nothing
        p = (K + 1) ÷ 2
    else
        @assert (K ÷ 2) ≤ p ≤ (K + 1) "The value of p must be between (K÷2) and (K+1)."
    end
    q = K + 1 - p
    return [hk[i + j - 1] for i in 1:p, j in 1:q]
end

    