"""
    takagi_factor(A::Matrix{ComplexF64}) -> (sv, Uz, err)

Perform the Takagi factorization on a complex symmetric matrix A.
Returns the vector of singular values `sv`, the matrix `Uz`,
and the norm of the difference between A and Uz*S*Uz^T.
"""

function takagi_factor(A::AbstractMatrix{T}; err::Bool=false) where T<:ComplexF64
    
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


function prony_sub(hk::Vector{ComplexF64}, eps::Float64)

    H = hankel_matrix(hk)
    
    # Takagi factorization
    sv, U = takagi_factor(H)
    # Determine the model order M
    M = count(>(eps * sv[1]), sv[1:end]) + 1
    
    # Roots
    uv = conj(U[:, M])
    roots = AMRVW.roots(uv)
    roots = sort(roots, by = x -> abs(x))
    γ = roots[1:M]

    return γ
end

function prony_sub(hk::Vector{ComplexF64}, M::Int)

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

function prony(hk::AbstractVector{<:ComplexF64}, dt::Real, eps::Real)
    gamm = prony_sub(hk, eps; ncols=ncols)
    expon, coeff = solve_vandermonde(hk, gamm, dt)
    return ExponentialFitting(expon, coeff)
end

function prony(func::Function, tmin::Float64, tmax::Float64, eps::Float64; nsamples::Int=501)
    # K must be odd
    @assert isodd(nsamples) "nsamples must be odd"
    hk = Vector{ComplexF64}(undef, nsamples)
    dt = (tmax-tmin) / (nsamples-1)
    hk = [func(dt * (k-1) + tmin) for k in 1:nsamples]
    return prony(hk, dt, eps; ncols=ncols)
end

function prony(hk::AbstractVector{<:ComplexF64}, dt::Real, M::Int)
    gamm = prony_sub(hk, eps; ncols=ncols)
    expon, coeff = solve_vandermonde(hk, gamm, dt)
    return ExponentialFitting(expon, coeff)
end

function prony(func::Function, tmin::Float64, tmax::Float64, M::Int; nsamples::Int=501)
    # K must be odd
    @assert isodd(nsamples) "nsamples must be odd for Prony method"
    hk = Vector{ComplexF64}(undef, nsamples)
    dt = (tmax-tmin) / (nsamples-1)
    hk = [func(dt * (k-1) + tmin) for k in 1:nsamples]
    return prony(hk, dt, M; ncols=ncols)
end


