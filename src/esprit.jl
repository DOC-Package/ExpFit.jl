"""
    esprit_sub(hk, eps; p=nothing)

Estimate the eigenvalues γ using the ESPRIT subspace method from the Hankel matrix
constructed from `hk`.
"""
function esprit_sub(hk::AbstractVector{<:ComplexF64}, eps::Real; cols::Union{Int,Nothing}=nothing)
    # Hankel matrix construction.
    H = hankel_matrix(hk; q=cols)

    # Perform SVD on the Hankel matrix.
    svd_res = svd(H, full=true)
    sv = svd_res.S
    V  = svd_res.V'
    
    # Determine the model order M as the number of singular values greater than eps * (largest singular value)
    M = count(>(eps * sv[1]), sv)
    
    cols = size(H, 2)
    W0 = transpose(V[1:M, 1:(cols-1)])
    W1 = transpose(V[1:M, 2:cols]) 
    FM = pinv(W0) * W1
    γ = eigvals(FM)

    return γ
end


"""
    esprit(hk, dt, eps; cols=nothing) -> exponent, coeff

Perform the ESPRIT algorithm using discrete data `hk` and the sampling interval `dt`.
"""
function esprit(hk::AbstractVector{<:ComplexF64}, dt::Real, eps::Real; cols::Union{Int,Nothing}=nothing)
    gamm = esprit_sub(hk, eps; cols=cols)
    # The function solve_vandermonde is assumed to be implemented in another file.
    expon, coeff = solve_vandermonde(hk, gamm, dt)
    return ExponentialFitting(expon, coeff)
end

"""
    esprit(func, tmin, tmax, K, eps; cols=nothing) -> (exponent, coeff)

Sample the function `func` over the interval [tmin, tmax] to create discrete data,
and then perform the ESPRIT algorithm.
"""
function esprit(func::Function, tmin::Real, tmax::Real, K::Int, eps::Real; cols::Union{Int,Nothing}=nothing)
    dt = (tmax - tmin) / (K - 1)
    hk = Vector{ComplexF64}(undef, K)
    hk = [func(tmin + dt * (k-1)) for k in 1:K]
    return esprit(hk, dt, eps; cols=cols)
end