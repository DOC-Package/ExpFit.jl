using ExpFit
using Test

@testset "ExpFit.jl" begin
    include("test_prony.jl")
    include("test_matrixpencil.jl")
    include("test_esprit.jl")
    include("test_espira1.jl")
    include("test_espira2.jl")
    include("test_fastesprit.jl")
    include("test_expfit.jl")
    include("test_baltru.jl")
end
