include("ftest.jl")
using Test
using LinearAlgebra
using ExpFit

@testset "balanced_truncation.jl" begin 

    tmin = 0.0
    tmax  = 50.0      
    eps = 1e-2
    N = 100
    t = range(tmin, tmax, length=N)
    
    a, c = generate_exponent_coefficient_pairs(100)
    f = Exponentials(a,c)

    ef = expred(a, c, eps)
    err = abs.(ef.(t) .- f.(t))
    @test norm(err)/sqrt(N) < eps

    ef = expred(a, c, 8)
    err = abs.(ef.(t) .- f.(t))
    @test norm(err)/sqrt(N) < eps

    a, c = generate_exponent_coefficient_pairs_real(100)
    f = Exponentials(a,c)

    ef = expred(a, c, eps)
    err = abs.(ef.(t) .- f.(t))
    @test norm(err)/sqrt(N) < eps

    ef = expred(a, c, 4)
    err = abs.(ef.(t) .- f.(t))
    @test norm(err)/sqrt(N) < eps

end
