include("../plot.jl")
using LinearAlgebra
using ExpFit
using Random

function generate_exponent_coefficient_pairs(n::Int)
    Random.seed!(123)
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

tmin = 0.0
tmax  = 50.0      
eps = 1e-2
N = 100
t = range(tmin, tmax, length=N)

a, c = generate_exponent_coefficient_pairs(100)
f = Exponentials(a,c)
er = expred(a, c, eps)
print("Approximation order = ", length(er.coeff), "\n")

t = range(tmin, tmax, length=N*2)
fv = f.(t)
erv = er.(t)
err = abs.(erv .- fv)
println("Root mean square = ", norm(err)/sqrt(N))
plot_res(t, fv, erv, err)