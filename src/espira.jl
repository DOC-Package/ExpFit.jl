"""
    espira1

This function first computes a modified FFT of y and places knots on the unit circle.
It then calls the AAA routine from RationalFunctionApproximation.jl to obtain a
rational approximant of (modified) DFT data, and finally extracts the ESPIRA parameters via a partial-fraction procedure.
"""
function espira1(f::Vector{ComplexF64}, dt::Real, eps::Real)

    N = length(f)
    # FFT
    F = fft(f)      

    # AAA
    Z = exp.(2π*im .* (0:N-1) ./ N)
    F = F .* (Z .^ (-1))
    bary = aaa(Z, F; tol=eps)
    
    # Extract the ESPIRA parameters
    pol = poles(bary)
    res = residues(bary, pol)
    coef = res ./ (1 .- pol.^N)
    freq = -log.(pol) ./ dt

    return Exponentials(freq, coef)
end

function espira1(func::Function, tmin::Real, tmax::Real, nsamples::Int, eps::Real)
    dt = (tmax - tmin) / (nsamples - 1)
    f = [func(tmin + dt * (k-1)) for k in 1:nsamples]
    return espira1(f, dt, eps)
end

function espira1(f::Vector{ComplexF64}, dt::Real, M::Int)

    N = length(f)
    # FFT
    F = fft(f)      

    # AAA
    Z = exp.(2π*im .* (0:M-1) ./ N)
    F = F .* (Z .^ (-1))
    bary = aaa(Z, F; max_degree=M)
    
    # Extract the ESPIRA parameters
    pol = poles(bary)
    res = residues(bary, pol)
    coef = res ./ (1 .- pol.^N)
    freq = -log.(pol) ./ dt

    return Exponentials(freq, coef)
end

function espira1(func::Function, tmin::Real, tmax::Real, nsamples::Int, M::Int)
    dt = (tmax - tmin) / (nsamples - 1)
    f = [func(tmin + dt * (k-1)) for k in 1:nsamples]
    return espira1(f, dt, M)
end

function espira2_sub!(z, f, Z, F, F1)

    # Compute f1 by matching each support point with its F1 value
    # and removing the support points from the approximant.
    f1 = similar(f)
    for (i, z_val) in enumerate(z)
        idx = findfirst(x -> isapprox(x, z_val; atol=1e-12), Z)
        f1[i] = F1[idx]
        deleteat!(Z, idx)
        deleteat!(F, idx)
        deleteat!(F1, idx)
    end

    # Form the Cauchy matrix
    m = length(z)
    C = hcat([1.0 ./ (Z .- z[j]) for j in 1:m]...)
    
    # Create the joint matrix L
    L0 = Diagonal(F) * C - C * Diagonal(f)
    L1 = Diagonal(F1) * C - C * Diagonal(f1)
    L = hcat(L0, L1)

    return L
end

"""
    espira2
    
This function computes a shifted FFT of its input data and uses AAA to obtain a rational approximant. 
It then extracts the exponential parameters by converting the resulting poles and residues into a sum of exponentials.
"""
function espira2(f::Vector{ComplexF64}, dt::Real, eps::Real)
    
    N = length(f)
    f0 = copy(f)
    # Compute the FFT of f.
    F = fft(f)     

    # AAA 
    Z = exp.(2π * im * (0:N-1) ./ N)
    F1 = copy(F)           
    F = F .* (Z .^ (-1))
    r = aaa(Z, F; tol=eps)

    # Create the joint matrix L
    L = espira2_sub!(r.nodes, r.values, Z, F, F1)

    # Perform SVD on the joint matrix L
    res = svd(L, full=true)
    sv = res.S
    W  = res.V'  
    
    # Determine the model order M 
    M = count(>(eps * sv[1]), sv) 
    m = length(r.nodes)
    W1 = W[1:M, 1:m]
    W2 = W[1:M, m+1:2m]
    γ = eigvals(pinv(transpose(W1)) * transpose(W2))

    expon, coeff = solve_vandermonde(f0, γ, dt)

    return Exponentials(expon, coeff)
end

function espira2(func::Function, tmin::Real, tmax::Real, nsamples::Int, eps::Real)
    dt = (tmax - tmin) / (nsamples - 1)
    f = [func(tmin + dt * (k-1)) for k in 1:nsamples]
    return espira2(f, dt, eps)
end

function espira2(f::Vector{ComplexF64}, dt::Real, M::Int)
    
    N = length(f)
    f0 = copy(f)
    # Compute the FFT of f.
    F = fft(f)     

    # AAA 
    Z = exp.(2π * im * (0:N-1) ./ N)
    F1 = copy(F)           
    F = F .* (Z .^ (-1))
    r = aaa(Z, F; max_degree=M)

    # Create the joint matrix L
    L = espira2_sub!(r.nodes, r.values, Z, F, F1)

    # Perform SVD on the joint matrix L
    res = svd(L, full=true)
    sv = res.S
    W  = res.V'  
    W1 = W[1:M, 1:m]
    W2 = W[1:M, m+1:2m]
    γ = eigvals(pinv(transpose(W1)) * transpose(W2))

    expon, coeff = solve_vandermonde(f0, γ, dt)

    return Exponentials(expon, coeff)
end

function espira2(func::Function, tmin::Real, tmax::Real, nsamples::Int, M::Int)
    dt = (tmax - tmin) / (nsamples - 1)
    f = [func(tmin + dt * (k-1)) for k in 1:nsamples]
    return espira2(f, dt, M)
end