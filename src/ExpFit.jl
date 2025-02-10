module ExpFit

export Exponentials, ESPRIT, MPencil, Prony
export prony, prony2, matrix_pencil, esprit, espira1, espira2, fast_esprit
export expfit, expred
export coneig

include("types.jl")
include("common.jl")
include("prony.jl")
include("matrix_pencil.jl")
include("esprit.jl")
include("espira.jl")
include("bidiagonalization.jl")
include("fast_esprit.jl")
include("pdcauchy.jl")
include("balanced_truncation.jl")

using LinearAlgebra
using AMRVW
using FFTW
using RationalFunctionApproximation
using Random

end
