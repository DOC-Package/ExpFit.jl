include("../plot.jl")
using Test
using LinearAlgebra
using ExpFit
using SpecialFunctions

tmin = 0.0
tmax  = 50.0      
eps = 1e-4     
N = 100    
t = range(tmin, tmax, length=N)
f = t -> besselj(0,t) + 1.0im*besselj(1,t)
ef = esprit(f, tmin, tmax, N, eps)
print("Approximation order = ", length(ef.coeff), "\n")

t = range(tmin, tmax, length=N*2)
fv = f.(t)
efv = ef.(t)
err = abs.(efv .- fv)
plot_res(t, fv, efv, err)
