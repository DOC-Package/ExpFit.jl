"""
    takagi_factor(A::Matrix{ComplexF64}) -> (sv, Uz, err)

Perform the Takagi factorization on a complex symmetric matrix A.
Returns the vector of singular values `sv`, the matrix `Uz`,
and the norm of the difference between A and Uz*S*Uz^T.
"""

function takagi_factor(A::AbstractMatrix{T}; err::Bool=false) where T<:Number
    
    n = size(A, 1)
    # check A is symmetric
    @assert isapprox(A, transpose(A), atol=1e-12) "A matrix must be symmetric for Takagi factorization"

    # Singular value decomposition
    svdres = svd(A, full=true)
    U  = svdres.U
    sv = collect(svdres.S)  
    V  = svdres.V  

    # construct Uz
    conjV = conj.(V) 
    Z = adjoint(U) * conjV   
    Zsr = sqrt(Z)
    Uz = U * Zsr

    # error
    if err
        A1 = Uz * Diagonal(sv) * transpose(Uz)
        err = norm(A1-A)
        return sv, Uz, err
    else
        return sv, Uz
    end
end


function prony_sub(hk::AbstractVector{<:Number}, eps::Float64)

    H = hankel_matrix(hk)
    
    # Takagi factorization
    sv, U = takagi_factor(H)
    # Determine the model order M
    M = count(>(eps * sv[1]), sv) 
    
    # Roots
    uv = conj(U[:, M])
    roots = AMRVW.roots(uv)
    roots = sort(roots, by = x -> abs(x))
    γ = roots[1:M]

    return γ
end

function prony_sub(hk::AbstractVector{<:Number}, M::Int)

    H = hankel_matrix(hk)
    
    # Takagi factorization
    sv, U = takagi_factor(H)
    
    # Roots
    uv = conj(U[:, M])
    roots = AMRVW.roots(uv)
    roots = sort(roots, by = x -> abs(x))
    γ = roots[1:M]

    return γ
end

"""
    prony(func, tmax, Nin, Mout)

Estimate the exponents and coefficients of the Prony series for the function `func`.    
"""

function prony(hk::AbstractVector{<:Number}, dt::Real, eps::Real)
    gamm = prony_sub(hk, eps)
    expon, coeff = solve_vandermonde(hk, gamm, dt)
    return Exponentials(expon, coeff)
end

function prony(func::Function, tmin::Float64, tmax::Float64, nsamples::Int, eps::Float64)
    @assert isodd(nsamples) "nsamples must be odd"
    dt = (tmax-tmin) / (nsamples-1)
    hk = [func(dt * (k-1) + tmin) for k in 1:nsamples]
    return prony(hk, dt, eps)
end

function prony(func::Function, tmin::Float64, tmax::Float64, dt::Real, eps::Float64)
    @assert isapprox((tmax-tmin)/dt, round((tmax-tmin)/dt), atol=1e-12) "(tmax-tmin)/dt must be an integer"
    nsamples = Int(round((tmax-tmin)/dt)) + 1
    hk = [func(dt * (k-1) + tmin) for k in 1:nsamples]
    return prony(hk, dt, eps)
end

function prony(hk::AbstractVector{<:Number}, dt::Real, M::Int)
    gamm = prony_sub(hk, M)
    expon, coeff = solve_vandermonde(hk, gamm, dt)
    return Exponentials(expon, coeff)
end

function prony(func::Function, tmin::Float64, tmax::Float64, nsamples::Int, M::Int)
    @assert isodd(nsamples) "nsamples must be odd for Prony method"
    dt = (tmax-tmin) / (nsamples-1)
    hk = [func(dt * (k-1) + tmin) for k in 1:nsamples]
    return prony(hk, dt, M)
end

function prony(func::Function, tmin::Float64, tmax::Float64, dt::Real, M::Int)
    @assert isapprox((tmax-tmin)/dt, round((tmax-tmin)/dt), atol=1e-12) "(tmax-tmin)/dt must be an integer"
    nsamples = Int(round((tmax-tmin)/dt)) + 1
    hk = [func(dt * (k-1) + tmin) for k in 1:nsamples]
    return prony(hk, dt, M)
end


function prony2(func::Function, tmin::Real, tmax::Real, nsamples::Int, eps::Real; ncols::Union{Int,Nothing}=nothing)
    
    dt = (tmax - tmin) / (nsamples - 1)
    hk = [func(tmin + dt * (k-1)) for k in 1:nsamples]

    H = hankel_matrix(hk[1:nsamples-1]; q=ncols)  # Note that this is the shifted Hankel matrix
    nrows, ncols = size(H)
    println("H ", size(H))

    hm = [hk[m + ncols - 1] for m in 1:nrows]
    q = pinv(H) * hm
    
    push!(q, 1.0)
    println("q ", size(q))
    roots = AMRVW.roots(q)
    println("roots ", size(roots))
    roots = sort(roots, by = x -> abs(x), rev=true)
    println("roots ", abs.(roots))
    γ = roots[1:ncols]

    coeff = solve_vandermonde!(hk, γ)

    M = ncols - count(x -> abs(x) < eps, coeff)
    println("M = $M")
    println("coeff", abs.(coeff))

    γ = γ[1:M]
    exponent, coeff = solve_vandermonde(hk, γ, dt)

    return Exponentials(exponent, coeff)
end