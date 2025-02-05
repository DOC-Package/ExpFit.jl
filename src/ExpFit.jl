module ExpFit

export ExponentialFitting
export prony, matrix_pencil, esprit

include("types.jl")
include("common.jl")
include("prony.jl")
include("matrix_pencil.jl")
include("esprit.jl")

using LinearAlgebra
using AMRVW

end
