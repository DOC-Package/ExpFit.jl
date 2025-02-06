module ExpFit

export ExponentialFitting, ESPRIT, MPencil, Prony
export prony, matrix_pencil, esprit, expfit
export coneig

include("types.jl")
include("common.jl")
include("prony.jl")
include("matrix_pencil.jl")
include("esprit.jl")

using LinearAlgebra
using AMRVW

end
