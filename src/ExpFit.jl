module ExpFit

export Exponentials, ESPRIT, MPencil, Prony
export prony, matrix_pencil, esprit, expfit, expred
export coneig

include("types.jl")
include("common.jl")
include("prony.jl")
include("matrix_pencil.jl")
include("esprit.jl")
include("pdcauchy.jl")
include("balanced_truncation.jl")

using LinearAlgebra
using AMRVW
#using FFTW
#using RationalFunctionApproximation

end
