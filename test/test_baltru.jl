using Test
using LinearAlgebra
using ExpFit
using SpecialFunctions
using Random

function generate_exponent_coefficient_pairs(n::Int)
    exponents = Vector{ComplexF64}(undef, n)
    coefficients = Vector{ComplexF64}(undef, n)
    for i in 1:n
        re_exp = rand() * 9.9 + 0.1
        im_exp = randn()  
        exponents[i] = re_exp + im_exp*im

        re_coeff = randn()
        im_coeff = randn()
        coefficients[i] = re_coeff + im_coeff*im
    end
    return exponents, coefficients
end

@testset "balanced_truncation.jl" begin 

    tmin = 0.0
    tmax  = 50.0      
    eps = 1e-2
    N = 100
    q = 33        
    
    # Create exponents and coefficients
    Random.seed!(1234)
    a, c = generate_exponent_coefficient_pairs(100)

    t = range(tmin, tmax, length=N)
    f = Exponentials(a,c)

    @time ef = expred(a, c, eps)
    err = abs.(ef.(t) .- f.(t))
    @test norm(err)/sqrt(N) < eps
    println("error: ", norm(err)/sqrt(N))
    println("degree: ", length(ef.expon))

    @time ef = expred(a, c, 8)
    err = abs.(ef.(t) .- f.(t))
    @test norm(err)/sqrt(N) < eps
    println("error: ", norm(err)/sqrt(N))
    println("degree: ", length(ef.expon))

end
