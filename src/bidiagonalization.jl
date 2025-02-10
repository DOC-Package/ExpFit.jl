function fast_fft_hankel_mv(L::Int, K::Int, fft_h::AbstractVector{<:Number}, x::AbstractVector{<:Number})
    @assert length(x) == K 

    # The minimum power of 2 greater than or equal to L+K-1
    P = nextpow(2, L + K - 1)

    # Extending the vector x
    tilde_x = vcat( reverse(x), zeros(eltype(x), P - K) )

    # Convolution calculation using FFT
    fft_x = fft(tilde_x)
    tilde_y = ifft( fft_h .* fft_x )

    return tilde_y[1:L]
end

function fast_fft_hankelH_mv(L::Int, K::Int, fft_h::AbstractVector{<:Number}, x::AbstractVector{<:Number})
    @assert length(x) == L 

    # The minimum power of 2 greater than or equal to L+K-1
    P = nextpow(2, L + K - 1)

    # Extending the vector x
    tilde_x = vcat( zeros(eltype(x), P - L), reverse(conj(x)) )

    # Convolution calculation using FFT
    fft_x = fft(tilde_x)
    tilde_y = ifft( fft_h .* fft_x )

    return tilde_y[P-K+1 : P]
end

function fft_hankel(L::Int, K::Int, h::AbstractVector{<:Number}) :: Vector{ComplexF64}
    @assert length(h) == L + K - 1

    # The minimum power of 2 greater than or equal to L+K-1
    P = nextpow(2, L + K - 1)

    # Extending the vector h
    tilde_h = vcat( h[K : L+K-1],
                    zeros(eltype(h), P - (L + K - 1)),
                    h[1:K-1] )

    # Convolution calculation using FFT
    fft_h = fft(tilde_h)

    return fft_h
end

"""
    fnrm(s)

Generic weighted norm computation
"""
function fnrm(s::AbstractVector{<:Number})
    n = length(s)
    n2 = div(n, 2)
    res1 = sum(i * abs2(s[i]) for i in 1:(n2 + 1))
    res2 = sum((n - i + 1) * abs2(s[i]) for i in (n2 + 2):n)
    return sqrt(res1 + res2)
end

"""
Generic Partial Lanczos Bidiagonalization

This function extends partial_lanczos_bidiagonalization_dble to handle both
real and complex input vectors s.
"""
function partial_lanczos_bidiagonalization(h::AbstractVector{T}, tol::Real) where T<:Number
    n = length(h)
    m = div(n, 2) + 1
    Fh = fft_hankel(m, m, h)
    FhH = fft_hankel(m, m, conj.(h))
    # Initialize qj
    qj = rand(m)
    qj[1] = 1.0 / sqrt(eps(Float64))
    qj ./= norm(qj, 2)
    Q = reshape(qj, m, 1)  
    # Compute initial left Lanczos vector pj
    pj = fast_fft_hankel_mv(m, m, Fh, qj)
    alpha = [norm(pj)]
    pj .= pj / alpha[1]
    P = reshape(pj, m, 1)  # Pmat holds left Lanczos vectors
    
    # Choose norm function based on type of s.
    nrm = fnrm(h)
    fnorm_A = alpha[1]^2
    epsj2 = 1.0
    j = 1
    
    println("Partial Lanczos Bidiagonalization")
    beta = []
    # Lanczos recursion loop
    while sqrt(abs(epsj2)) / nrm > tol
        # Recursion for right Lanczos vector:
        qj1 = fast_fft_hankelH_mv(m, m, FhH, pj) - alpha[j] * qj
        qj1 .-= Q * (Q' * qj1)  # Orthogonalize against previous q's
        b_j = norm(qj1)
        push!(beta, b_j)
        qj1 ./= b_j
        Q = hcat(Q, qj1)
        
        # Recursion for left Lanczos vector:
        pj1 = fast_fft_hankel_mv(m, m, Fh, qj1) - b_j * pj
        pj1 .-= P * (P' * pj1)  # Orthogonalize against previous p's
        a_j1 = norm(pj1, 2)
        push!(alpha, a_j1)
        pj1 ./= a_j1
        P = hcat(P, pj1)
        
        # Update for next iteration:
        qj = copy(qj1)
        pj = copy(pj1)
        t = a_j1^2 + b_j^2
        fnorm_A += t
        epsj2 = t / fnorm_A
        j += 1
        println("Iteration $j: residual = ", sqrt(abs(epsj2)) / nrm)
    end
    return j, alpha, beta, P, Q
end

function partial_lanczos_bidiagonalization(h::AbstractVector{T}, step::Int) where T<:Number
    n = length(h)
    m = div(n, 2) + 1
    Fh = fft_hankel(m, m, h)
    FhH = fft_hankel(m, m, conj.(h))
    # Initialize qj
    qj = rand(m)
    qj[1] = 1.0 / sqrt(eps(Float64))
    qj ./= norm(qj, 2)
    Q = reshape(qj, m, 1)  
    # Compute initial left Lanczos vector pj
    pj = fast_fft_hankel_mv(m, m, Fh, qj)
    alpha = [norm(pj)]
    pj .= pj / alpha[1]
    P = reshape(pj, m, 1)  # Pmat holds left Lanczos vectors
    beta = []
    # Lanczos recursion loop
    for j = 1:step
        # Recursion for right Lanczos vector:
        qj1 = fast_fft_hankelH_mv(m, m, FhH, pj) - alpha[j] * qj
        qj1 .-= Q * (Q' * qj1)  # Orthogonalize against previous q's
        b_j = norm(qj1)
        push!(beta, b_j)
        qj1 ./= b_j
        Q = hcat(Q, qj1)
        
        # Recursion for left Lanczos vector:
        pj1 = fast_fft_hankel_mv(m, m, Fh, qj1) - b_j * pj
        pj1 .-= P * (P' * pj1)  # Orthogonalize against previous p's
        a_j1 = norm(pj1, 2)
        push!(alpha, a_j1)
        pj1 ./= a_j1
        P = hcat(P, pj1)
        
        # Update for next iteration:
        qj = copy(qj1)
        pj = copy(pj1)
    end
    return alpha, beta, P, Q
end
