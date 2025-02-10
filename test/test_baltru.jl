include("ftest.jl")
using Test
using LinearAlgebra
using ExpFit

@testset "balanced_truncation.jl" begin 

    tmin = 0.0
    tmax  = 50.0      
    eps = 1e-2
    N = 100
    
    # Create exponents and coefficients
    a, c = generate_exponent_coefficient_pairs(100)

    t = range(tmin, tmax, length=N)
    f = Exponentials(a,c)

    @time ef = expred(a, c, eps)
    err = abs.(ef.(t) .- f.(t))
    @test norm(err)/sqrt(N) < eps
    println("error: ", norm(err)/sqrt(N))
    println("degree: ", length(ef.expon))

    @time ef = expred(a, c, 8)
    err = abs.(ef.(t) .- f.(t))
    @test norm(err)/sqrt(N) < eps
    println("error: ", norm(err)/sqrt(N))
    println("degree: ", length(ef.expon))

end
