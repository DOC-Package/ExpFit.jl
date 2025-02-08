module ExpFit

export Exponentials, ESPRIT, MPencil, Prony
export prony, matrix_pencil, esprit, espira1, espira2
export expfit, expred
export coneig

include("types.jl")
include("common.jl")
include("prony.jl")
include("matrix_pencil.jl")
include("esprit.jl")
include("espira.jl")
include("pdcauchy.jl")
include("balanced_truncation.jl")

using LinearAlgebra
using AMRVW
using FFTW
using RationalFunctionApproximation

end
