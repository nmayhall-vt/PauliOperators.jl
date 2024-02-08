using PauliOperators
using Test

@testset "PauliOperators.jl" begin
    include("tests.jl")
    include("test_matvec.jl")
    include("test_vectorized.jl")
end
