function esprit_sub(hk::AbstractVector{<:Number}, eps::Real; ncols::Union{Int,Nothing}=nothing)
    # Hankel matrix construction.
    H = hankel_matrix(hk; q=ncols)

    # Perform SVD on the Hankel matrix.
    svd_res = svd(H, full=true)
    sv = svd_res.S
    V  = svd_res.V'
    
    # Determine the model order M 
    M = count(>(eps * sv[1]), sv)
    
    q = size(H, 2)
    W0 = transpose(V[1:M, 1:(q-1)])
    W1 = transpose(V[1:M, 2:q]) 
    FM = pinv(W0) * W1
    γ = eigvals(FM)

    return γ
end

function esprit_sub(hk::AbstractVector{<:Number}, M::Int; ncols::Union{Int,Nothing}=nothing)
    # Hankel matrix construction.
    H = hankel_matrix(hk; q=ncols)

    # Perform SVD on the Hankel matrix.
    svd_res = svd(H, full=true)
    V  = svd_res.V'
    
    q = size(H, 2)
    W0 = transpose( V[1:M, 1:(q-1)] )
    W1 = transpose( V[1:M, 2:q] ) 
    FM = pinv(W0) * W1
    γ = eigvals(FM)

    return γ
end

"""
    esprit(hk::AbstractVector{<:Number}, dt::Real, eps::Real; ncols::Union{Int,Nothing}=nothing) :: Exponentials

Perform the ESPRIT algorithm using discrete data `hk` and the sampling interval `dt` for a given tolerance `eps`.
"""
function esprit(hk::AbstractVector{<:Number}, dt::Real, eps::Real; ncols::Union{Int,Nothing}=nothing)
    gamm = esprit_sub(hk, eps; ncols=ncols)
    expon, coeff = solve_vandermonde(hk, gamm, dt)
    return Exponentials(expon, coeff)
end

"""
    esprit(func::Function, tmin::Real, tmax::Real, nsamples::Int, eps::Real; ncols::Union{Int,Nothing}=nothing) :: Exponentials

Perform the ESPRIT algorithm using a function `func` in the range [tmin,tmax] and `nsamples` sampling points for a given tolerance `eps`.
"""
function esprit(func::Function, tmin::Real, tmax::Real, nsamples::Int, eps::Real; ncols::Union{Int,Nothing}=nothing)
    dt = (tmax - tmin) / (nsamples - 1)
    hk = [func(tmin + dt * (k-1)) for k in 1:nsamples]
    return esprit(hk, dt, eps; ncols=ncols)
end

"""
    esprit(func::Function, tmin::Real, tmax::Real, dt::Real, eps::Real; ncols::Union{Int,Nothing}=nothing) :: Exponentials

Perform the ESPRIT algorithm using a function `func` in the range [tmin,tmax] and a sampling interval `dt` for a given tolerance `eps`.
"""
function esprit(func::Function, tmin::Real, tmax::Real, dt::Real, eps::Real; ncols::Union{Int,Nothing}=nothing)
    @assert isapprox((tmax-tmin)/dt, round((tmax-tmin)/dt), atol=1e-12) "(tmax-tmin)/dt must be an integer"
    nsamples = Int(round((tmax-tmin)/dt)) + 1
    hk = [func(tmin + dt * (k-1)) for k in 1:nsamples]
    return esprit(hk, dt, eps; ncols=ncols)
end

"""
    esprit(hk::Vector{<:Number}, dt::Real, M::Int; ncols::Union{Int,Nothing}=nothing) :: Exponentials

Perform the ESPRIT algorithm using discrete data `hk` and the sampling interval `dt` for a given model order `M`.
"""
function esprit(hk::AbstractVector{<:Number}, dt::Real, M::Int; ncols::Union{Int,Nothing}=nothing)
    gamm = esprit_sub(hk, M; ncols=ncols)
    expon, coeff = solve_vandermonde(hk, gamm, dt)
    return Exponentials(expon, coeff)
end

"""
    esprit(func::Function, tmin::Real, tmax::Real, nsamples::Int, M::Int; ncols::Union{Int,Nothing}=nothing) :: Exponentials

Perform the ESPRIT algorithm using a function `func` in the range [tmin,tmax] and `nsamples` sampling points for a given model order `M`.
"""
function esprit(func::Function, tmin::Real, tmax::Real, nsamples::Int, M::Int; ncols::Union{Int,Nothing}=nothing)
    dt = (tmax - tmin) / (nsamples - 1)
    hk = [func(tmin + dt * (k-1)) for k in 1:nsamples]
    return esprit(hk, dt, M; ncols=ncols)
end

"""
    esprit(func::Function, tmin::Real, tmax::Real, dt::Real, M::Int; ncols::Union{Int,Nothing}=nothing) :: Exponentials

Perform the ESPRIT algorithm using a function `func` in the range [tmin,tmax] and a sampling interval `dt` for a given model order `M`.
"""
function esprit(func::Function, tmin::Real, tmax::Real, dt::Real, M::Int; ncols::Union{Int,Nothing}=nothing)
    @assert isapprox((tmax-tmin)/dt, round((tmax-tmin)/dt), atol=1e-12) "(tmax-tmin)/dt must be an integer"
    nsamples = Int(round((tmax-tmin)/dt)) + 1
    hk = [func(tmin + dt * (k-1)) for k in 1:nsamples]
    return esprit(hk, dt, M; ncols=ncols)
end