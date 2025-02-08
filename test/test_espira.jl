using Test
using LinearAlgebra
using ExpFit
using SpecialFunctions

@testset "espira.jl" begin 

    tmin = 0.0
    tmax  = 50.0      
    eps = 1e-4    
    N = 100
    q = 33             

    t = range(tmin, tmax, length=N)
    f = t -> besselj(0,t) + 1.0im*besselj(1,t)

    @time ef = espira2(f, tmin, tmax, N, eps)
    println("degree: ", length(ef.coeff))
    err = abs.(ef.(t) .- f.(t))
    @test norm(err)/sqrt(N) < eps
    println("error: ", norm(err)/sqrt(N))

end
