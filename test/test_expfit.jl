include("ftest.jl")
using Test
using LinearAlgebra
using ExpFit

@testset "expfit" begin 

    tmin = 0.0
    tmax  = 20.0      
    eps = 1e-4    
    N = 101             

    t = range(tmin, tmax, length=N)
    dt = t[2] - t[1]

    a, c = generate_exponent_coefficient_pairs(100)
    f = Exponentials(a,c)
    fv = f.(t)
    fmax = maximum(abs.(fv))

    ef = expfit(f, tmin, tmax, N, eps)
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax

    ef = expfit(f, tmin, tmax, N, eps; alg=Pencil())
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax

    ef = expfit(f, tmin, tmax, N, eps; alg=Prony())
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax*2

    ef = expfit(f, tmin, tmax, N, eps; alg=ESPIRA1())
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax

    ef = expfit(f, tmin, tmax, N, eps; alg=ESPIRA2())
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax

    a, c = generate_exponent_coefficient_pairs_real(100)
    f = Exponentials(a,c)
    fv = f.(t)
    fmax = maximum(abs.(fv))

    ef = expfit(f, tmin, tmax, N, eps; alg=FastESPRIT())
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax

end
