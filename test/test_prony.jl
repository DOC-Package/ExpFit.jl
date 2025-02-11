include("ftest.jl")
using Test
using LinearAlgebra
using ExpFit

@testset "prony.jl" begin 

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

    ef = prony(f, tmin, tmax, N, eps)
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax*10

    ef = prony(f, tmin, tmax, N, 12)
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax*10

    ef = prony(f, tmin, tmax, dt, eps)
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax*10

    ef = prony(f, tmin, tmax, dt, 12)
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax*10

    ef = prony(fv, dt, eps)
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax*10

    ef = prony(fv, dt, 12)
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax*10

    a, c = generate_exponent_coefficient_pairs_real(100)
    f = Exponentials(a,c)
    fv = f.(t)
    fmax = maximum(abs.(fv))

    ef = prony(f, tmin, tmax, N, eps)
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax*10

    ef = prony(f, tmin, tmax, N, 12)
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax*10

    ef = prony(f, tmin, tmax, dt, eps)
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax*10

    ef = prony(f, tmin, tmax, dt, 12)
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax*10

    ef = prony(fv, dt, eps)
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax*10

    ef = prony(fv, dt, 12)
    err = abs.(ef.(t) .- fv)
    @test norm(err)/sqrt(N) < eps*fmax*10

end
