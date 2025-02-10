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
    M = count(>(eps * sv[1]), sv) + 1
    
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
    if iseven(length(hk))
        hk = hk[1:end-1]
        println("nsamples must be odd for Prony method, removing the last element")
    end
    gamm = prony_sub(hk, eps)
    expon, coeff = solve_vandermonde(hk, gamm, dt)
    return Exponentials(expon, coeff)
end

function prony(func::Function, tmin::Float64, tmax::Float64, nsamples::Int, eps::Float64)
    if iseven(nsamples)
        nsamples += 1
        println("nsamples must be odd for Prony method, adding 1")
    end
    dt = (tmax-tmin) / (nsamples-1)
    hk = [func(dt * (k-1) + tmin) for k in 1:nsamples]
    return prony(hk, dt, eps)
end

function prony(func::Function, tmin::Float64, tmax::Float64, dt::Real, eps::Float64)
    @assert isapprox((tmax-tmin)/dt, round((tmax-tmin)/dt), atol=1e-12) "(tmax-tmin)/dt must be an integer"
    nsamples = Int(round((tmax-tmin)/dt)) + 1
    if iseven(nsamples)
        nsamples += 1
        println("nsamples must be odd for Prony method, adding 1")
    end
    hk = [func(dt * (k-1) + tmin) for k in 1:nsamples]
    return prony(hk, dt, eps)
end

function prony(hk::AbstractVector{<:Number}, dt::Real, M::Int)
    if iseven(length(hk))
        hk = hk[1:end-1]
        println("nsamples must be odd for Prony method, removing the last element")
    end
    gamm = prony_sub(hk, M)
    expon, coeff = solve_vandermonde(hk, gamm, dt)
    return Exponentials(expon, coeff)
end

function prony(func::Function, tmin::Float64, tmax::Float64, nsamples::Int, M::Int)
    if iseven(nsamples)
        nsamples += 1
        println("nsamples must be odd for Prony method, adding 1")
    end
    dt = (tmax-tmin) / (nsamples-1)
    hk = [func(dt * (k-1) + tmin) for k in 1:nsamples]
    return prony(hk, dt, M)
end

function prony(func::Function, tmin::Float64, tmax::Float64, dt::Real, M::Int)
    @assert isapprox((tmax-tmin)/dt, round((tmax-tmin)/dt), atol=1e-12) "(tmax-tmin)/dt must be an integer"
    nsamples = Int(round((tmax-tmin)/dt)) + 1
    if iseven(nsamples)
        nsamples += 1
        println("nsamples must be odd for Prony method, adding 1")
    end
    hk = [func(dt * (k-1) + tmin) for k in 1:nsamples]
    return prony(hk, dt, M)
end