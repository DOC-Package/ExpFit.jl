"""
    matrix_pencil_sub(hk, L, eps) -> gamm

Estimate the eigenvalues γ using the Matrix Pencil method from the Hankel matrix
constructed from the discrete data `hk`.

- `hk`    : A vector of complex samples. It is assumed that the length of `hk` is 2N.
- `L`     : An integer parameter for constructing the Hankel matrix. The Hankel matrix
            is built with dimensions (2N - L) × (L + 1).
- `eps` : A threshold used for model order determination. Specifically, the model order
            M is determined as the smallest index m (with 1 ≤ m ≤ L) such that 
            |R[m+1, m+1]|/|R[1,1]| < eps. In that case, M is set to m+1;
            if no such index is found, M is set to 1.

Returns: A tuple `gamm` where:
- `gamm` is a vector of eigenvalues obtained from the pencil.
- `M` is the detected model order.
"""
function matrix_pencil_sub(hk::AbstractVector{<:ComplexF64}, eps::Real; q::Union{Int,Nothing}=nothing)
    # Construct the Hankel matrix.
    H = hankel_matrix(hk; q=q)
    
    # Perform QR decomposition with column pivoting.
    #res = qr(H, Val(true))
    res = qr(H, ColumnNorm())
    R = Matrix(res.R) * transpose(res.P)
    
    
    # Determine the model order M by inspecting the diagonal of R.
    # Find the first m (1 ≤ m ≤ L) for which |R[m+1, m+1]|/|R[1,1]| < eps.
    L = size(H, 2) - 1
    m_opt = findfirst(m -> abs(R[m + 1, m + 1]) / abs(R[1, 1]) < eps, 1:L)
    M = m_opt !== nothing ? m_opt + 1 : 1

    # Form the sub-blocks (transposed) of R to construct the pencil.
    T0t = transpose(R[1:M, 1:L])
    T1t = transpose(R[1:M, 2:(L+1)])
    FM = pinv(T0t) * T1t
    gamm = eigvals(FM)

    return gamm
end


"""
    matrix_pencil(hk, tmin, tmax, L, eps) -> (exponent, coeff)

Perform the Matrix Pencil method using discrete data `hk` and the sampling interval defined by [tmin, tmax].

- `hk`    : Discrete data (a vector of complex numbers), expected to be of length 2N.
- `tmin`  : The start of the sampling interval.
- `tmax`  : The end of the sampling interval.
- `L`     : Parameter for the Hankel matrix construction (the matrix will be of size (2N-L)×(L+1)).
- `eps` : Threshold for determining the model order in the Matrix Pencil method.

Returns: A tuple `(exponent, coeff)` where:
- `exponent` is a vector of estimated exponents (e.g. growth/decay rates).
- `coeff`    is a vector of estimated coefficients.
  
Note: The Vandermonde system is solved via least squares.
"""

"""
function matrix_pencil(hk::AbstractVector{<:ComplexF64}, dt::Real, eps::Real; q::Union{Int,Nothing}=nothing)
    gamm = matrix_pencil_sub(hk, eps; q=q)
    # The function solve_vandermonde is assumed to be implemented in another file.
    exponent, coeff = solve_vandermonde(hk, gamm, dt)
    
    return exponent, coeff
end
"""

"""
    matrix_pencil(func, tmin, tmax, K, L, eps) -> (exponent, coeff)

Sample the function `func` over the interval [tmin, tmax] at K equally spaced points,
and then perform the Matrix Pencil method.

- `func`  : The function to be sampled. It should accept a real argument `t` and return a ComplexF64.
- `tmin`  : The start of the sampling interval.
- `tmax`  : The end of the sampling interval.
- `K`     : Number of samples to generate (must be 2N, i.e. even).
- `L`     : Parameter for the Hankel matrix construction.
- `eps` : Threshold for determining the model order.

Returns: A tuple `(exponent, coeff)` containing the estimated exponents and coefficients.
"""

"""
function matrix_pencil(func::Function, tmin::Real, tmax::Real, K::Int, eps::Real; q::Union{Int,Nothing}=nothing)
    hk = Vector{ComplexF64}(undef, 2N)
    dt = (tmax - tmin) / (K - 1)
    hk = [func(tmin + dt * (k - 1)) for k in 1:K]
    return matrix_pencil(hk, dt, eps; q=q)
end
"""

function matrix_pencil(
    func::Function,
    tmin::Real,
    tmax::Real,  
    N::Int,  
    eps::Real;       
    q::Union{Int,Nothing}=nothing
)

    # Generate discrete samples
    hk = Vector{ComplexF64}(undef, N)
    dt = (tmax - tmin) / (N - 1)
    hk = [func(tmin + dt * k) for k in 0:N-1]

    """
    # Hankel matrix
    #Hkm = Matrix{ComplexF64}(undef, 2N-L, L+1)
    #Hkm = [hk[k + m - 1] for k in 1:(2N - L), m in 1:(L + 1)]

    Hkm = hankel_matrix(hk; q=L+1)

    # QR decomposition
    res = qr(Hkm, ColumnNorm())
    R = Matrix(res.R) * transpose(res.P)

    # Determine M
    m_opt = findfirst(m -> abs(R[m + 1, m + 1]) / abs(R[1, 1]) < epsin, 1:L)
    M = m_opt !== nothing ? m_opt + 1 : 1

    T0t = transpose(R[1:M, 1:L])
    T1t = transpose(R[1:M, 2:(L+1)])
    FM = pinv(T0t) * T1t
    gamm = eigvals(FM)
    """

    gamm = matrix_pencil_sub(hk, eps; q=q)

    exponent, coeff = solve_vandermonde(hk, gamm, dt)

    return exponent, coeff
end


function matrix_pencil(
    func::Function,
    tmin::Real,
    tmax::Real,  
    eps::Real,  
    N::Int,          
    L::Int   
)

    # Generate discrete samples
    hk = Vector{ComplexF64}(undef, 2N)
    dt = (tmax - tmin) / (2N - 1)
    hk = [func(tmin + dt * k) for k in 0:2N-1]

    gamm = matrix_pencil_sub(hk, eps; q=L+1)

    exponent, coeff = solve_vandermonde(hk, gamm, dt)

    return exponent, coeff
end