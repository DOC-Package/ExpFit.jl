##############################################################################
# espira.jl
#
# Julia version of the ESPIRA method.
#
# This code computes an exponential sum approximation from a given vector y.
#
# It uses the RationalFunctionApproximation.jl package to perform AAA.
#
##############################################################################

using FFTW                    # for fft
using LinearAlgebra           # for eigen, cond, etc.
using RationalFunctionApproximation  # provides the AAA algorithm

"""
    espira1appr(y, Nexp)

Input:
  • y    : vector of samples.
  • Nexp : desired order (number of exponentials) in the constructed sum.

Output:
  • g    : vector of coefficients in the exponential sum.
  • pol  : vector of poles (from the partial-fraction representation)
           of the rational function.
           
This function first computes a modified FFT of y and places knots on the unit circle.
It then calls the AAA routine from RationalFunctionApproximation.jl to obtain a
rational approximant of (modified) DFT data, and finally extracts the ESPIRA parameters
via a partial-fraction procedure (implemented in `prz` below).
"""
function espira1appr(y, Nexp)
    # Compute the FFT of y.
    F = fft(y)
    M = length(F)          # number of modified DFT values

    # Set knots on the unit circle.
    # (Here the knots are Z = exp(2πi * k/M), for k = 0,1,...,M-1.)
    Z = exp.(2π*im .* (0:M-1) ./ M)
    
    # Make F and Z into vectors (if not already) and modify F:
    # F = F .* Z.^(-1)
    F = F .* (Z .^ (-1))
    
    # ---- Use AAA from RationalFunctionApproximation.jl -----------------------
    # We now approximate the function defined on the knots Z with data F.
    # (In the MATLAB code the AAA iteration is run for m = 2:Nexp+1 support points.
    # Here we request a rational approximant of degree Nexp.)
    approx = aaa(F, Z; degree = Nexp)
    
    # The AAA result (here “approx”) is assumed to have fields:
    #   approx.support         :: vector of support points (denoted z in MATLAB)
    #   approx.f               :: corresponding data values at support points (f)
    #   approx.w               :: barycentric weights (w)
    #
    # (Adjust the field names if necessary.)
    z = approx.support
    fvals = approx.f
    w = approx.w
    
    # Remove support points with (nearly) zero weight.
    tolw = 1e-13
    keep = findall(x -> abs(x) > tolw, w)
    z = z[keep]
    fvals = fvals[keep]
    w = w[keep]
    
    # ---- Compute the ESPIRA parameters via a partial-fraction procedure ----
    # In the MATLAB code this is done by the function "prz". The following Julia
    # function computes:
    #   • pol   : the poles of the rational approximant (via a generalized eigenvalue problem)
    #   • AB    : intermediate partial-fraction coefficients,
    #   • g     : the final exponential sum coefficients,
    # with g = –AB ./ (1 – pol.^M).
    g, AB, pol = prz(F, Z, z, fvals, w, M)
    
    return g, pol
end

##############################################################################
# prz
#
# Computes the partial-fraction decomposition of the rational approximant.
#
# Input:
#   • F : the modified DFT values (unused in the internal calculation but passed along)
#   • Z : the vector of knots on the unit circle.
#   • z : the support points (from AAA)
#   • f : the corresponding data values at the support points.
#   • w : the barycentric weights from AAA.
#   • M : the original length (used to compute the correction term 1 – pol^M).
#
# Output:
#   • g   : vector of coefficients in the exponential sum.
#   • AB  : vector of intermediate partial-fraction coefficients.
#   • pol : vector of poles computed via a generalized eigenvalue problem.
##############################################################################
function prz(F, Z, z, f, w, M)
    m = length(w)
    # Build (m+1)x(m+1) matrices B and E.
    # B is the identity with the (1,1) entry replaced by 0.
    B = Matrix{ComplexF64}(I, m+1, m+1)
    B[1,1] = 0
    
    # E is constructed so that its first row is [0, w₁, w₂, …, wₘ]
    # and rows 2:end have 1 in the first column and the diagonal entry zᵢ in the (i+1,i+1) slot.
    E = zeros(ComplexF64, m+1, m+1)
    for j in 1:m
        E[1, j+1] = w[j]
    end
    for i in 2:(m+1)
        E[i, 1] = 1
        E[i, i] = z[i-1]
    end
    
    # Solve the generalized eigenvalue problem E * v = λ B * v.
    # (The MATLAB code then takes eigenvalues that are finite.)
    ge = eigen(E, B)
    pol_all = ge.values
    pol = [λ for λ in pol_all if isfinite(λ)]
    
    # Build the Cauchy matrix CC with entries: CC(i,j) = -1/(z[i] - pol[j]).
    m2 = length(z)
    nPol = length(pol)
    CC = Array{ComplexF64}(undef, m2, nPol)
    for i in 1:m2
        for j in 1:nPol
            CC[i,j] = -1 / (z[i] - pol[j])
        end
    end
    
    # Solve for the partial-fraction coefficients AB from CC * AB = f.
    AB = CC \ f
    
    # Compute the exponential sum coefficients:
    #   g = -AB ./ (1 - pol.^M)
    g = similar(AB)
    for j in 1:nPol
        g[j] = -AB[j] / (1 - pol[j]^M)
    end
    
    # (For diagnostic purposes, print the condition number of CC.)
    cond_CC = cond(CC)
    println("Condition number of the Cauchy matrix for ESPIRA I: ", cond_CC)
    
    return g, AB, pol
end

# ===================== Example usage =====================
# (Uncomment the block below to test the ESPIRA approximation.)

# using Random
# Random.seed!(1234)
#
# # Create a sample signal y (for example, a noisy sum of exponentials)
# M = 128
# t = 0:M-1
# y = exp.(0.1im .* t) .+ 0.5*exp.(-0.2im .* t) .+ 0.1*randn(M)
#
# # Choose the desired exponential sum order, e.g. Nexp = 5.
# Nexp = 5
#
# # Compute the ESPIRA approximation.
# g, pol = espira1appr(y, Nexp)
#
# println("Exponential sum coefficients g:")
# println(g)
# println("Poles:")
# println(pol)
