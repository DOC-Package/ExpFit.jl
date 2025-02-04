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

function calc_roots(a::Vector{Complex{T}}) where T<:Real
    nd = length(a)
    
    Ca = zeros(Complex{T}, nd, nd)
    for n in 1:(nd-1)
        Ca[n+1, n] = 1.0 + 0im
    end

    for n in 1:nd
        Ca[n, nd] = -a[n]
    end

    # Roots
    roots = eigen(Ca).values

    # check the error
    czero = Vector{Complex{T}}(undef, nd)
    for n in 1:nd
        czero[n] = roots[n]^nd
        for m in 1:nd
            czero[n] += a[m] * roots[n]^(m-1)
        end
    end
    err = norm(czero)

    return roots, err
end


"""
    prony(func, tmax, Nin, Mout)

# Arguments
- `func::Function`
    - 実数 t を受け取り複素数を返す関数 (Fortran の `complex(real64), external :: func`)
- `tmax::Float64`
    - サンプリング区間 [0, tmax]
- `Nin::Int`
    - サンプリング・行列次元などに使う整数パラメータ
- `Mout::Int`
    - 出力モード数 (Fortran で intent(inout) だったもの)

# Returns
- `(Mout_new, rhoout, alphaout)`
"""
function prony(
    func::Function,
    tmin::Float64,
    tmax::Float64,
    eps::Float64,
    N::Int
)

    hk = Vector{ComplexF64}(undef, 2N+1)
    dt = (tmax-tmin) / (2N)
    hk = [func(dt * k + tmin) for k in 0:2N]

    Hkl = Matrix{ComplexF64}(undef, N+1, N+1)
    Hkl = [hk[k+l-1] for k in 1:N+1, l in 1:N+1]
    
    # Takagi factorization
    sv, U, err = takagi_factor(Hkl)
    M = count(>(eps * sv[1]), sv[1:N+1]) + 1

    uvt = copy(conj(U[:, M]))  
    uvt ./= uvt[N+1]      
    uv = copy(uvt[1:N])
    roots, err = calc_roots(uv) 
    println("err: ", err)
    roots = sort(roots, by = x -> abs(x))
    gamm = roots[1:M]

    # Solve the Vandermonde system
    Gmn = Matrix{ComplexF64}(undef, 2N+1, M)
    Gmn = [gamm[m]^k for k in 0:2N, m in 1:M]

    # Coefficients
    coeff = Gmn \ hk

    # Exponents
    exponent = Vector{ComplexF64}(undef, M)
    exponent = [-log(gamm[m]) / dt for m in 1:M]

    # Sort by magnitude
    idx = sortperm(coeff, by = x -> abs(x), rev=true)
    coeff = coeff[idx]
    exponent = gamm[idx]

    return exponent, coeff
end
