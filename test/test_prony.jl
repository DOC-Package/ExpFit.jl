include("utils.jl")
using Test
using LinearAlgebra
using ExpFit
using SpecialFunctions

@testset "prony.jl" begin 

    tmin = 0.0
    tmax  = 50.0      
    eps = 1e-3     
    N = 101
    q = 33             

    t = range(tmin, tmax, length=2*N)
    f = t -> besselj(0,t) + 1.0im*besselj(1,t)

    @time ef = prony(f, tmin, tmax, eps; nsamples=N)
    err = abs.(ef.(t) .- f.(t))
    @test norm(err) < eps*10.0
    print("error = ", norm(err), "\n")

end
