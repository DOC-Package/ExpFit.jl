"""
    matrix_pencil_sub(hk, eps; cols) -> γ

Estimate the eigenvalues γ using the Matrix Pencil method from the Hankel matrix
constructed from the discrete data `hk`.
"""
function matrix_pencil_sub(hk::AbstractVector{<:ComplexF64}, eps::Real; cols::Union{Int,Nothing}=nothing)
    # Construct the Hankel matrix.
    H = hankel_matrix(hk; q=cols)
    
    # Perform QR decomposition with column pivoting.
    res = qr(H, ColumnNorm())
    R = Matrix(res.R) * transpose(res.P)
    
    # Determine the model order M by inspecting the diagonal of R.
    # Find the first m (1 ≤ m ≤ L) for which |R[m+1, m+1]|/|R[1,1]| < eps.
    q = size(H, 2)
    m_opt = findfirst(m -> abs(R[m + 1, m + 1]) / abs(R[1, 1]) < eps, 1:(q-1)) + 1
    M = m_opt !== nothing ? m_opt + 1 : 1

    # Form the sub-blocks (transposed) of R to construct the pencil.
    T0t = transpose(R[1:M, 1:q-1])
    T1t = transpose(R[1:M, 2:q])
    FM = pinv(T0t) * T1t
    γ = eigvals(FM)

    return γ
end


"""
    matrix_pencil_sub(hk, dt, eps; cols) -> γ

Perform the Matrix Pencil method using discrete data `hk` and the sampling interval defined by [tmin, tmax].
"""
function matrix_pencil(hk::AbstractVector{<:ComplexF64}, dt::Real, eps::Real; cols::Union{Int,Nothing}=nothing)
    gamm = matrix_pencil_sub(hk, eps; cols=cols)
    # The function solve_vandermonde is assumed to be implemented in another file.
    exponent, coeff = solve_vandermonde(hk, gamm, dt)
    
    return exponent, coeff
end


"""
    matrix_pencil_sub(func, tmin, tmax, K, eps; cols) -> γ

Perform the Matrix Pencil method using discrete data `hk` and the sampling interval defined by [tmin, tmax].
"""
function matrix_pencil(func::Function, tmin::Real, tmax::Real, K::Int, eps::Real; cols::Union{Int,Nothing}=nothing)
    dt = (tmax - tmin) / (K - 1)
    hk = [func(tmin + dt * (k - 1)) for k in 1:K]
    return matrix_pencil(hk, dt, eps; cols=cols)
end
