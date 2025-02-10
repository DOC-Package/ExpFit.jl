"""
    solve_vandermonde()

Solve the overdetermined Vandermode system
"""

function solve_vandermonde(
    hk::Vector{<:Number},
    γ::Vector{<:Number},
    dt::Float64
)
    
    K = length(hk)
    M = length(γ)
    
    # Solve the Vandermonde system
    G = Matrix{ComplexF64}(undef, K, M)
    G = [γ[m]^k for k in 0:K-1, m in 1:M]

    # Coefficients
    coeff = G \ hk

    # Exponents
    exponent = Vector{ComplexF64}(undef, M)
    exponent .= -log.(γ) ./ dt
    
    # Sort by magnitude
    idx = sortperm(coeff, by = x -> abs(x), rev=true)
    exponent = exponent[idx]
    coeff = coeff[idx]

    return exponent, coeff
end

function solve_vandermonde!(
    hk::Vector{ComplexF64},
    γ::Vector{ComplexF64}
)
    
    K = length(hk)
    M = length(γ)
    
    # Solve the Vandermonde system
    G = Matrix{ComplexF64}(undef, K, M)
    G = [γ[m]^k for k in 0:K-1, m in 1:M]

    # Coefficients
    coeff = G \ hk
    
    # Sort by magnitude
    idx = sortperm(coeff, by = x -> abs(x), rev=true)
    γ = γ[idx]
    coeff = coeff[idx]

    return coeff
end

"""
    hankel_matrix(hk; q=nothing)

Construct a Hankel matrix from the discrete data `hk` (a vector of complex numbers).
"""
function hankel_matrix(hk::AbstractVector{<:Number}; q::Union{Int,Nothing}=nothing) :: Matrix{ComplexF64}
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