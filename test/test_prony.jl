using Test
using LinearAlgebra
using ExpFit
using SpecialFunctions

@testset "prony.jl" begin 

    tmin = 0.0
    tmax  = 50.0      
    eps = 1e-4     
    N = 201             

    t = range(tmin, tmax, length=N)
    dt = t[2] - t[1]
    f = t -> besselj(0,t) + 1.0im*besselj(1,t)

    """
    @time ef = prony(f, tmin, tmax, N, eps)
    err = abs.(ef.(t) .- f.(t))
    @test norm(err)/sqrt(N) < eps*2
    print("error = ", norm(err)/sqrt(N))

    @time ef = prony(f, tmin, tmax, dt, eps)
    err = abs.(ef.(t) .- f.(t))
    @test norm(err)/sqrt(N) < eps*2
    print("error = ", norm(err)/sqrt(N))
    """
    N = 800
    @time ef = prony2(f, tmin, tmax, N, eps)
    err = abs.(ef.(t) .- f.(t))
    @test norm(err)/sqrt(N) < eps*2
    print("error = ", norm(err)/sqrt(N))

end
