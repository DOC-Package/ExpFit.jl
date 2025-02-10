include("ftest.jl")
using Test
using LinearAlgebra
using ExpFit

@testset "fast_esprit.jl" begin 

    tmin = 0.0
    tmax  = 20.0      
    eps = 1e-6     
    N = 101           

    t = range(tmin, tmax, length=N)
    a, c = generate_exponent_coefficient_pairs_real(100)
    f = Exponentials(a,c)

    ef = fast_esprit(f, tmin, tmax, N, eps)
    err = abs.(ef.(t) .- f.(t))
    @test norm(err)/sqrt(N) < eps*10
    println("error: ", norm(err)/sqrt(N))

    """
    ef = expfit(f, tmin, tmax, N, eps; alg=ESPRIT())
    err = abs.(ef.(t) .- f.(t))
    @test norm(err)/sqrt(N) < eps
    #println("error: ", norm(err)/sqrt(N))

    
    @time exponent, coeff = esprit(f, tmin, tmax, N, eps; cols=q)
    err = [abs(sumexp(ti,exponent,coeff) - f(ti)) for ti in t]
    @test norm(err) < eps*10.0

    f = t -> besselj(2,t) + 1.0im*besselj(3,t)
    dt = t[2] - t[1]
    f_disc = [f(ti) for ti in t]
    @time exponent, coeff = esprit(f_disc, dt, eps)
    err = [abs(sumexp(ti,exponent,coeff) - f(ti)) for ti in t]
    @test norm(err) < eps*10.0

    tmin = 0.0
    tmax  = 50.0      
    eps = 1e-4     
    N = 105
    q = 35             

    t = range(tmin, tmax, length=N)
    f = t -> besselj(0,t) + 1.0im*besselj(1,t)

    @time exponent, coeff = esprit(f, tmin, tmax, N, eps)
    err = [abs(sumexp(ti,exponent,coeff) - f(ti)) for ti in t]
    @test norm(err) < eps*10.0

    @time exponent, coeff = esprit(f, tmin, tmax, N, eps; cols=q)
    err = [abs(sumexp(ti,exponent,coeff) - f(ti)) for ti in t]
    @test norm(err) < eps*10.0

    f = t -> besselj(2,t) + 1.0im*besselj(3,t)
    dt = t[2] - t[1]
    f_disc = [f(ti) for ti in t]
    @time exponent, coeff = esprit(f_disc, dt, eps)
    err = [abs(sumexp(ti,exponent,coeff) - f(ti)) for ti in t]
    @test norm(err) < eps*10.0
    """
end
