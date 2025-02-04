include("utils.jl")
using Test
using LinearAlgebra
using ExpFit
using SpecialFunctions

@testset "matrix_pencil.jl" begin 

    tmin = 0.0
    tmax  = 50.0      
    eps = 1e-4     
    N = 100
    q = 33             

    t = range(tmin, tmax, length=2*N)
    f = t -> besselj(0,t) + 1.0im*besselj(1,t)

    @time exponent, coeff = matrix_pencil(f, tmin, tmax, N, eps)
    err = [abs(f_approx(ti,exponent,coeff) - f(ti)) for ti in t]
    @test norm(err) < eps*10.0

    @time exponent, coeff = matrix_pencil(f, tmin, tmax, N, eps; q=q)
    err = [abs(f_approx(ti,exponent,coeff) - f(ti)) for ti in t]
    @test norm(err) < eps*10.0

    f = t -> besselj(2,t) + 1.0im*besselj(3,t)
    dt = t[2] - t[1]
    f_disc = [f(ti) for ti in t]
    @time exponent, coeff = matrix_pencil(f_disc, dt, eps)
    err = [abs(f_approx(ti,exponent,coeff) - f(ti)) for ti in t]
    @test norm(err) < eps*10.0

    tmin = 0.0
    tmax  = 50.0      
    eps = 1e-4     
    N = 105
    q = 35             

    t = range(tmin, tmax, length=N)
    f = t -> besselj(0,t) + 1.0im*besselj(1,t)

    @time exponent, coeff = matrix_pencil(f, tmin, tmax, N, eps)
    err = [abs(f_approx(ti,exponent,coeff) - f(ti)) for ti in t]
    @test norm(err) < eps*10.0

    @time exponent, coeff = matrix_pencil(f, tmin, tmax, N, eps; q=q)
    err = [abs(f_approx(ti,exponent,coeff) - f(ti)) for ti in t]
    @test norm(err) < eps*10.0

    f = t -> besselj(2,t) + 1.0im*besselj(3,t)
    dt = t[2] - t[1]
    f_disc = [f(ti) for ti in t]
    @time exponent, coeff = matrix_pencil(f_disc, dt, eps)
    err = [abs(f_approx(ti,exponent,coeff) - f(ti)) for ti in t]
    @test norm(err) < eps*10.0
end
