include("utils.jl")
using Test
using LinearAlgebra
using ExpFit
using SpecialFunctions

@testset "prony.jl" begin 

    tmin = 0.0
    tmax  = 50.0      
    eps = 1e-3     
    N = 50
    q = 33             

    t = range(tmin, tmax, length=2*N)
    f = t -> besselj(0,t) + 1.0im*besselj(1,t)

    @time exponent, coeff = prony(f, tmin, tmax, N, eps)
    err = [abs(sumexp(ti,exponent,coeff) - f(ti)) for ti in t]
    print("error =", norm(err), "\n")
    @test norm(err) < eps*10.0

end
