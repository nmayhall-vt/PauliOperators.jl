using PauliOperators
using Test

@testset "Paulis.jl" begin
    include("test_operator_methods.jl")
    include("test_Pauli.jl")
    include("test_Ket.jl")
end
