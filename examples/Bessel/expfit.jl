include("../plot.jl")
using Test
using LinearAlgebra
using ExpFit
using SpecialFunctions

tmin = 0.0
tmax  = 50.0      
eps = 1e-3     
N = 100    
t = range(tmin, tmax, length=N)
f = t -> besselj(0,t) + 1.0im*besselj(1,t)
ef = expfit(f, tmin, tmax, N, eps)
print("Approximation order = ", length(ef.coeff), "\n")

t = range(tmin, tmax, length=N*2)
fv = f.(t)
efv = ef.(t)
err = abs.(efv .- fv)
println("Root mean square = ", norm(err)/sqrt(N))
plot_res(t, fv, efv, err)
