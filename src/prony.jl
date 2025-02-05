"""
    takagi_factor(A::Matrix{ComplexF64}) -> (sv, Uz, err)

Perform the Takagi factorization on a complex symmetric matrix A.
Returns the vector of singular values `sv`, the matrix `Uz`,
and the norm of the difference between A and Uz*S*Uz^T.
"""

function takagi_factor(A::AbstractMatrix)
    
    n = size(A, 1)
    # check A is symmetric
    @assert isapprox(A, transpose(A), atol=1e-12) "A must be symmetric"

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
    S = Diagonal(sv)
    A1 = Uz * S * transpose(Uz)
    err = norm(A1-A)

    return sv, Uz, err
end

"""
    calc_roots(a::Vector{Complex{T}}) where T<:Real

Calculate the roots of a polynomial with coefficients `a`.
"""

function calc_roots(a::Vector{ComplexF64}) 
    nd = length(a)-1
    a = a[1:nd] ./ a[end]
    
    Ca = zeros(ComplexF64, nd, nd)
    for n in 1:(nd-1)
        Ca[n+1, n] = 1.0 + 0im
    end

    for n in 1:nd
        Ca[n, nd] = -a[n]
    end

    # Roots
    roots = eigen(Ca).values

    # check the error
    czero = Vector{ComplexF64}(undef, nd)
    for n in 1:nd
        czero[n] = roots[n]^nd
        for m in 1:nd
            czero[n] += a[m] * roots[n]^(m-1)
        end
    end
    err = norm(czero)

    return roots, err
end

function roots_error(coeff::Vector{ComplexF64}, roots::Vector{ComplexF64})
    nd = length(roots)
    czero = Vector{ComplexF64}(undef, nd)
    for n in 1:nd
        czero[n] = roots[n]^nd
        for m in 1:nd
            czero[n] += coeff[m] * roots[n]^(m-1)
        end
    end
    err = norm(czero)

    return err
end


"""
    prony(func, tmax, Nin, Mout)

Estimate the exponents and coefficients of the Prony series for the function `func`.    
"""
function prony(
    func::Function,
    tmin::Float64,
    tmax::Float64,
    N::Int,
    eps::Float64
)

    hk = Vector{ComplexF64}(undef, 2N+1)
    dt = (tmax-tmin) / (2N)
    hk = [func(dt * k + tmin) for k in 0:2N]

    Hkl = Matrix{ComplexF64}(undef, N+1, N+1)
    Hkl = [hk[k+l-1] for k in 1:N+1, l in 1:N+1]
    
    # Takagi factorization
    sv, U, err = takagi_factor(Hkl)
    @assert err < 1e-12 "error in Takagi factorization: $err"
    
    M = count(>(eps * sv[1]), sv[1:N+1]) + 1

    uv = conj(U[:, M])
    
    # Roots
    roots = AMRVW.roots(uv)
    roots = sort(roots, by = x -> abs(x))
    γ = roots[1:M]

    expon, coeff = solve_vandermonde(hk, γ, dt)

    return expon, coeff
end
