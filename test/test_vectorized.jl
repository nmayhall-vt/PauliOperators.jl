using PauliOperators
using LinearAlgebra
using Random
using Test
using BlockDavidson
# using LinearMaps
using BenchmarkTools

@testset "vectorized" begin
   
    N = 10

    a = Pauli(N)
    H = a + a'

    for i in 1:40
        ai = rand()*random_Pauli(N)
        H += ai + ai'
    end


    ρ = random_Pauli(N) + random_Pauli(N)
   
    ρv = VectorizedPauliSum(ρ)
    Hr = VectorizedRMult(H)
    Hl = VectorizedLMult(H)
    HU = VectorizedConjugate(H)
    HC = VectorizedCommutator(H)

    ρmat = Matrix(ρ)
    Hmat = Matrix(H)

    a = PauliOperators.commutator(ρ,H)
    display(a)

    aref = ρmat * Hmat - Hmat * ρmat

    err = abs(norm(Matrix(a) - aref))
    @show err     
    @test abs(err) < 1e-14     

    err = norm(Matrix(VectorizedCommutator(H)*VectorizedPauliSum(ρ)) + aref)
    @test abs(err) < 1e-14     

    # Need to finish testing the functions in type_vectorized...
end
# run()