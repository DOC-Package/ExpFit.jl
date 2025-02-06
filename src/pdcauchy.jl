function mytril(A::AbstractMatrix{T}) where T<:Number
    m, n = size(A)
    L = zeros(ComplexF64, m, n)
    for i in 2:m
        jmax = min(i - 1, n)
        L[i, 1:jmax] .= A[i, 1:jmax]
    end
    return L
end

function pivot_order!(a::AbstractVector{T}, b::AbstractVector{T}, x::AbstractVector{T}, y::AbstractVector{T}) where T<:Number
    n = length(a)
    g = [ a[i]*b[i] / (x[i] + y[i]) for i in 1:n ]

    # initialize permutation matrix as identity
    P = Matrix{T}(I, n, n)
    for m in 1:(n-1)
        l = argmax(abs.(@view g[m:end])) + m - 1
        g[m], g[l] = g[l], g[m]
        x[m], x[l] = x[l], x[m]
        y[m], y[l] = y[l], y[m]
        a[m], a[l] = a[l], a[m]
        b[m], b[l] = b[l], b[m]
        P[:, m], P[:, l] = P[:, l], P[:, m]
        @views g[m+1:end] .*= ((x[m+1:end] .- x[m]) ./ (x[m+1:end] .+ y[m])) .*
                              ((y[m+1:end] .- y[m]) ./ (y[m+1:end] .+ x[m]))
    end
    return P
end

function partial_cholesky!(a::AbstractVector{T}, b::AbstractVector{T}, x::AbstractVector{T}, y::AbstractVector{T}) where T<:Number
    n = length(a)
    P = pivot_order!(a, b, x, y)
    G = zeros(T, n, n)
    @views G[:, 1] .= a .* b[1] ./ (x .+ y[1])
    for k in 2:n
        @views begin
            a[k:end] .*= (x[k:end] .- x[k-1]) ./ (x[k:end] .+ y[k-1])
            b[k:end] .*= (y[k:end] .- y[k-1]) ./ (y[k:end] .+ x[k-1])
            G[k:end, k] .= a[k:end] .* b[k] ./ (x[k:end] .+ y[k])
        end
    end
    D = Diagonal(sqrt.(diag(G))) |> Matrix
    D2i = Diagonal(1.0 ./ diag(G)) |> Matrix
    L = mytril(G) * D2i
    @views L[diagind(L)] .+= one(T)

    return L, D, P
end

function coneig_rrd(X::Matrix{T}, D::Matrix{T}) where T<:Number
    m = size(X, 1)
    G = D * (transpose(X) * (X * D))
    
    # QR factorization of G
    F = qr(G, ColumnNorm())
    R = Matrix(F.R) * transpose(Matrix(F.P))
    
    # SVD of R
    SVD = svd(R); s = SVD.S; Us = SVD.U
    
    invD = Diagonal(1 ./ diag(D))
    R1 = invD * R * invD
    X1 = invD * Us * Diagonal(sqrt.(s)) 
    Y1 = R1 \ X1
    U = conj(X * Y1)
    
    # phase correction
    for i in 1:m
        ti = sum(U[:, i].^2)
        phase = ti / abs(ti)
        U[:, i] .*= sqrt(conj(phase))
    end

    return s, U
end

function coneig(a::Vector{T}, b::Vector{T}, x::Vector{T}, y::Vector{T}) where T<:Number
    L, D, P = partial_cholesky!(map(copy, (a, b, x, y))...)
    Xt = P * L
    s, U = coneig_rrd(Xt, D)
    return s, U
end

