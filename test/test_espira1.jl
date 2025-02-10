include("ftest.jl")
using Test
using LinearAlgebra
using ExpFit

@testset "espira1.jl" begin 

    tmin = 0.0
    tmax  = 20.0      
    eps = 1e-4    
    N = 100             

    t = range(tmin, tmax, length=N)
    dt = t[2] - t[1]

    a, c = generate_exponent_coefficient_pairs(100)
    f = Exponentials(a,c)
    fv = f.(t)
    fmax = maximum(abs.(fv))

    ef = espira1(f, tmin, tmax, N, eps)
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax

    ef = espira1(f, tmin, tmax, N, 9)
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax

    ef = espira1(f, tmin, tmax, dt, eps)
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax

    ef = espira1(f, tmin, tmax, dt, 9)
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax

    ef = espira1(fv, dt, eps)
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax

    ef = espira1(fv, dt, 9)
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax

    a, c = generate_exponent_coefficient_pairs_real(100)
    f = Exponentials(a,c)
    fv = f.(t)
    fmax = maximum(abs.(fv))

    ef = espira1(f, tmin, tmax, N, eps)
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax

    ef = espira1(f, tmin, tmax, N, 9)
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax

    ef = espira1(f, tmin, tmax, dt, eps)
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax

    ef = espira1(f, tmin, tmax, dt, 9)
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax

    ef = espira1(fv, dt, eps)
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax

    ef = espira1(fv, dt, 9)
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax

end
