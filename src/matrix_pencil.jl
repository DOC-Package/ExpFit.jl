function matrix_pencil_sub(hk::AbstractVector{<:Number}, eps::Real; ncols::Union{Int,Nothing}=nothing)
    # Construct the Hankel matrix.
    H = hankel_matrix(hk; q=ncols)
    
    # Perform QR decomposition with column pivoting.
    res = qr(H, ColumnNorm())
    S = Matrix(res.R) * transpose(res.P)
    
    # Determine the model order M by inspecting the diagonal of R.
    q = size(H, 2)
    M = findfirst(m -> abs(res.R[m + 1, m + 1]) / abs(res.R[1, 1]) < eps, 1:(q-1))
    
    # Form the sub-blocks (transposed) of R to construct the pencil.
    T0t = transpose(S[1:M, 1:q-1])
    T1t = transpose(S[1:M, 2:q])
    FM = pinv(T0t) * T1t
    γ = eigvals(FM)

    return γ
end

function matrix_pencil_sub(hk::AbstractVector{<:Number}, M::Int; ncols::Union{Int,Nothing}=nothing)
    # Construct the Hankel matrix.
    H = hankel_matrix(hk; q=ncols)
    
    # Perform QR decomposition with column pivoting.
    res = qr(H, ColumnNorm())
    S = Matrix(res.R) * transpose(res.P)
    
    q = size(H, 2)
    # Form the sub-blocks (transposed) of R to construct the pencil.
    T0t = transpose(S[1:M, 1:q-1])
    T1t = transpose(S[1:M, 2:q])
    FM = pinv(T0t) * T1t
    γ = eigvals(FM)

    return γ
end

""" 
    matrix_pencil(hk::AbstractVector{<:Number}, dt::Real, eps::Real; ncols::Union{Int,Nothing}=nothing) :: Exponentials
    
Perform the matrix pencil algorithm using discrete data `hk` and the sampling interval `dt` for a given tolerance `eps`.
"""
function matrix_pencil(hk::AbstractVector{<:Number}, dt::Real, eps::Real; ncols::Union{Int,Nothing}=nothing)
    gamm = matrix_pencil_sub(hk, eps; ncols=ncols)
    exponent, coeff = solve_vandermonde(hk, gamm, dt)
    return Exponentials(exponent, coeff)
end

"""
    matrix_pencil(func::Function, tmin::Real, tmax::Real, nsamples::Int, eps::Real; ncols::Union{Int,Nothing}=nothing)

Perform the matrix pencil algorithm using a function `func` in the range [tmin,tmax] and `nsamples` sampling points for a given tolerance `eps`.
"""
function matrix_pencil(func::Function, tmin::Real, tmax::Real, nsamples::Int, eps::Real; ncols::Union{Int,Nothing}=nothing)
    dt = (tmax - tmin) / (nsamples - 1)
    hk = [func(tmin + dt * (k - 1)) for k in 1:nsamples]
    return matrix_pencil(hk, dt, eps; ncols=ncols)
end

"""
    matrix_pencil(func::Function, tmin::Real, tmax::Real, dt::Real, eps::Real; ncols::Union{Int,Nothing}=nothing) :: Exponentials

Perform the matrix pencil algorithm using a function `func` in the range [tmin,tmax] and a sampling interval `dt` for a given tolerance `eps`.
"""
function matrix_pencil(func::Function, tmin::Real, tmax::Real, dt::Real, eps::Real; ncols::Union{Int,Nothing}=nothing)
    @assert isapprox((tmax-tmin)/dt, round((tmax-tmin)/dt), atol=1e-12) "(tmax-tmin)/dt must be an integer"
    nsamples = Int(round((tmax-tmin)/dt)) + 1
    hk = [func(tmin + dt * (k - 1)) for k in 1:nsamples]
    return matrix_pencil(hk, dt, eps; ncols=ncols)
end

"""
    matrix_pencil(hk::Vector{<:Number}, dt::Real, M::Int; ncols::Union{Int,Nothing}=nothing) :: Exponentials

Perform the matrix pencil algorithm using discrete data `hk` and the sampling interval `dt` for a given model order `M`.
"""
function matrix_pencil(hk::AbstractVector{<:Number}, dt::Real, M::Int; ncols::Union{Int,Nothing}=nothing)
    gamm = matrix_pencil_sub(hk, M; ncols=ncols)
    exponent, coeff = solve_vandermonde(hk, gamm, dt)
    return Exponentials(exponent, coeff)
end

"""
    matrix_pencil(func::Function, tmin::Real, tmax::Real, nsamples::Int, M::Int; ncols::Union{Int,Nothing}=nothing) :: Exponentials

Perform the matrix pencil algorithm using a function `func` in the range [tmin,tmax] and `nsamples` sampling points for a given model order `M`.
"""
function matrix_pencil(func::Function, tmin::Real, tmax::Real, nsamples::Int, M::Int; ncols::Union{Int,Nothing}=nothing)
    dt = (tmax - tmin) / (nsamples - 1)
    hk = [func(tmin + dt * (k - 1)) for k in 1:nsamples]
    return matrix_pencil(hk, dt, M; ncols=ncols)
end

"""
    matrix_pencil(func::Function, tmin::Real, tmax::Real, dt::Real, M::Int; ncols::Union{Int,Nothing}=nothing) :: Exponentials

Perform the matrix pencil algorithm using a function `func` in the range [tmin,tmax] and a sampling interval `dt` for a given model order `M`.
"""
function matrix_pencil(func::Function, tmin::Real, tmax::Real, dt::Real, M::Int; ncols::Union{Int,Nothing}=nothing)
    @assert isapprox((tmax-tmin)/dt, round((tmax-tmin)/dt), atol=1e-12) "(tmax-tmin)/dt must be an integer"
    nsamples = Int(round((tmax-tmin)/dt)) + 1
    hk = [func(tmin + dt * (k - 1)) for k in 1:nsamples]
    return matrix_pencil(hk, dt, M; ncols=ncols)
end
