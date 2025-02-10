"""
    fast_esprit_sub(hk, eps; p=nothing)

Estimate the eigenvalues γ using the ESPRIT subspace method from the Hankel matrix
constructed from `hk`.
"""
function fast_esprit_sub(hk::AbstractVector{<:Number}, eps::Real; ncols::Union{Int,Nothing}=nothing)
    # Hankel matrix construction.
    H = hankel_matrix(hk; q=ncols)

    # Perform SVD on the Hankel matrix.
    M, alpha, beta, P, Q = partial_lanczos_bidiagonalization(hk, eps)
    
    q = size(H, 2)
    Q0 = Q[1:(q-1),1:M]
    Q1 = Q[2:q, 1:M]
    FM = pinv(Q0) * Q1
    γ = eigvals(FM)

    return γ
end

function fast_esprit_sub(hk::AbstractVector{<:Number}, M::Int; ncols::Union{Int,Nothing}=nothing)
    # Hankel matrix construction.
    H = hankel_matrix(hk; q=ncols)

    # Perform SVD on the Hankel matrix.
    alpha, beta, P, Q = partial_lanczos_bidiagonalization(hk, M)
    
    q = size(H, 2)
    Q0 = Q[1:M, 1:(q-1)]
    Q1 = Q[1:M, 2:q]
    FM = pinv(Q0) * Q1
    γ = eigvals(FM)

    return γ
end


"""
    fast_esprit(hk, dt, eps; ncols=nothing) -> exponent, coeff

Perform the ESPRIT algorithm using discrete data `hk` and the sampling interval `dt`.
"""
function fast_esprit(hk::AbstractVector{<:Number}, dt::Real, eps::Real; ncols::Union{Int,Nothing}=nothing)
    gamm = fast_esprit_sub(hk, eps; ncols=ncols)
    expon, coeff = solve_vandermonde(hk, gamm, dt)
    return Exponentials(expon, coeff)
end

function fast_esprit(func::Function, tmin::Real, tmax::Real, nsamples::Int, eps::Real; ncols::Union{Int,Nothing}=nothing)
    dt = (tmax - tmin) / (nsamples - 1)
    hk = [func(tmin + dt * (k-1)) for k in 1:nsamples]
    return fast_esprit(hk, dt, eps; ncols=ncols)
end

function fast_esprit(func::Function, tmin::Real, tmax::Real, dt::Real, eps::Real; ncols::Union{Int,Nothing}=nothing)
    @assert isapprox((tmax-tmin)/dt, round((tmax-tmin)/dt), atol=1e-12) "(tmax-tmin)/dt must be an integer"
    nsamples = Int(round((tmax-tmin)/dt)) + 1
    hk = [func(tmin + dt * (k-1)) for k in 1:nsamples]
    return fast_esprit(hk, dt, eps; ncols=ncols)
end

function fast_esprit(hk::AbstractVector{<:Number}, dt::Real, M::Int; ncols::Union{Int,Nothing}=nothing)
    gamm = fast_esprit_sub(hk, M; ncols=ncols)
    expon, coeff = solve_vandermonde(hk, gamm, dt)
    return Exponentials(expon, coeff)
end

function fast_esprit(func::Function, tmin::Real, tmax::Real, nsamples::Int, M::Int; ncols::Union{Int,Nothing}=nothing)
    dt = (tmax - tmin) / (nsamples - 1)
    hk = [func(tmin + dt * (k-1)) for k in 1:nsamples]
    return fast_esprit(hk, dt, M; ncols=ncols)
end

function fast_esprit(func::Function, tmin::Real, tmax::Real, dt::Real, M::Int; ncols::Union{Int,Nothing}=nothing)
    @assert isapprox((tmax-tmin)/dt, round((tmax-tmin)/dt), atol=1e-12) "(tmax-tmin)/dt must be an integer"
    nsamples = Int(round((tmax-tmin)/dt)) + 1
    hk = [func(tmin + dt * (k-1)) for k in 1:nsamples]
    return fast_esprit(hk, dt, M; ncols=ncols)
end