"""
    esprit_sub(hk, eps; p=nothing)

Estimate the eigenvalues γ using the ESPRIT subspace method from the Hankel matrix
constructed from `hk`.

- `hk` : A vector of complex data.
- `eps`: Threshold used to distinguish noise. Singular values smaller than `eps` times the largest singular value are considered noise.
- `p`  : Number of rows for the Hankel matrix (optional).

Returns: A vector of eigenvalues γ.
"""
function esprit_sub(hk::AbstractVector{<:ComplexF64}, eps::Real; p::Union{Int,Nothing}=nothing)
    H = hankel_matrix(hk; p=p)
    svd_res = svd(H, full=true)
    sv = svd_res.S
    V  = svd_res.V'
    
    # Determine the model order M as the number of singular values greater than eps * (largest singular value)
    M = count(>(eps * sv[1]), sv)
    
    q = size(H, 2)
    @assert q > 1 "Insufficient number of samples. The Hankel matrix must have more than one column."
    W0 = transpose( V[1:M, 1:(q-1)] )
    W1 = transpose( V[1:M, 2:q] ) 
    FM = pinv(W0) * W1
    
    gamm = eigvals(FM)
    return gamm
end


"""
    esprit(hk, dt, eps; p=nothing) -> (exponent, coeff)

Perform the ESPRIT algorithm using discrete data `hk` and the sampling interval `dt`.

- `hk` : Discrete data (vector of complex numbers).
- `dt` : Sampling interval.
- `eps`: Threshold for singular value determination.
- `p`  : Number of rows for the Hankel matrix (optional).

Returns: A tuple containing the estimated exponents and coefficients.
Note: The solution of the Vandermonde system is delegated to the function `solve_vandermonde`, which is assumed to be implemented elsewhere.
"""
function esprit(hk::AbstractVector{<:ComplexF64}, dt::Real, eps::Real; p::Union{Int,Nothing}=nothing)
    gamm = esprit_sub(hk, eps; p=p)
    # The function solve_vandermonde is assumed to be implemented in another file.
    exponent, coeff = solve_vandermonde(hk, gamm, dt)
    return exponent, coeff
end

"""
    esprit(func, tmin, tmax, N, eps; p=nothing) -> (exponent, coeff)

Sample the function `func` over the interval [tmin, tmax] to create discrete data,
and then perform the ESPRIT algorithm.

- `func` : The function to be sampled (should return a real or complex number).
- `tmin, tmax` : Endpoints of the sampling interval.
- `N`    : Number of samples (this implementation generates 2N data points).
- `eps`  : Threshold for singular value determination.
- `p`    : Number of rows for the Hankel matrix (optional).

Returns: A tuple containing the estimated exponents and coefficients.
"""
function esprit(func::Function, tmin::Real, tmax::Real, K::Int, eps::Real; p::Union{Int,Nothing}=nothing)
    #K = 2 * N  # Adjust to 2N-1 if needed.
    dt = (tmax - tmin) / (K - 1)
    hk = [func(tmin + dt * (k-1)) for k in 1:K]
    return esprit(hk, dt, eps; p=p)
end