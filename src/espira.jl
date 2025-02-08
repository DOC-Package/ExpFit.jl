"""
    espira1(y, Nexp)

This function first computes a modified FFT of y and places knots on the unit circle.
It then calls the AAA routine from RationalFunctionApproximation.jl to obtain a
rational approximant of (modified) DFT data, and finally extracts the ESPIRA parameters
via a partial-fraction procedure (implemented in `prz` below).
"""
function espira1(func::Function, tmin::Real, tmax::Real, nsamples::Int, eps::Real; ncols::Union{Int,Nothing}=nothing)
    dt = (tmax - tmin) / (nsamples - 1)
    f = [func(tmin + dt * (k-1)) for k in 1:nsamples]
    # Compute the FFT of f.
    F = fft(f)     

    Z = exp.(2π*im .* (0:nsamples-1) ./ nsamples)
    F = F .* (Z .^ (-1))
    
    bary = aaa(Z, F; tol=eps)
    
    # Extract the ESPIRA parameters
    pol = poles(bary)
    res = residues(bary, pol)
    M = length(F)
    coef = res ./ (1 .- pol.^M)
    freq = -log.(pol) ./ dt

    return Exponentials(freq, coef)
end

function espira1(f::Vector{ComplexF64}, eps::Real; ncols::Union{Int,Nothing}=nothing)
    dt = (tmax - tmin) / (nsamples - 1)
    f = [func(tmin + dt * (k-1)) for k in 1:nsamples]
    # Compute the FFT of y.
    F = fft(f)
    M = length(F)        

    Z = exp.(2π*im .* (0:M-1) ./ M)
    F = F .* (Z .^ (-1))
    
    bary = aaa(Z, F; tol=eps)
    
    # Extract the ESPIRA parameters
    pol = poles(bary)
    res = residues(bary, pol)
    coef = res ./ (1 .- pol.^M)
    freq = -log.(pol) ./ dt

    return Exponentials(freq, coef)
end


function espira2(func::Function, tmin::Real, tmax::Real, nsamples::Int, eps::Real; ncols::Union{Int,Nothing}=nothing)
    N = nsamples
    dt = (tmax - tmin) / (nsamples - 1)
    f = [func(tmin + dt * (k-1)) for k in 1:nsamples]
    f0 = copy(f)

    # Compute the FFT of f.
    F = fft(f)   
    Z = exp.(2π * im * (0:nsamples-1) ./ nsamples)
    F1 = copy(F)           
    F = F .* (Z .^ (-1))  

    # AAA routine from RationalFunctionApproximation.jl
    jmax = N ÷ 2
    r = aaa(Z, F; tol=eps)
    z = r.nodes
    f = r.values

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

    # Perform SVD on the joint matrix L
    res = svd(L, full=true)
    sv = res.S
    W  = res.V'  
    
    # Determine the model order M 
    M = count(>(eps * sv[1]), sv) 
    W1 = W[1:M, 1:m]
    W2 = W[1:M, m+1:2m]
    γ = eigvals(pinv(transpose(W1)) * transpose(W2))

    expon, coeff = solve_vandermonde(f0, γ, dt)

    return Exponentials(expon, coeff)
end
