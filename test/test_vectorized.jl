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

    ρmat = Matrix(ρ)
    Hmat = Matrix(H)

    a = PauliOperators.commutator(ρ,H)
    display(a)

    aref = ρmat * Hmat - Hmat * ρmat

    err = abs(norm(Matrix(a) - aref))
    @show err     
    @test abs(err) < 1e-14     

    # Commutator

    L = PauliOperators.vectorized_commutator(H)

    # display(L)

    @test abs(err) < 1e-14     
    a = L*ρv 

    err = norm(Matrix(a.ps) + aref)
    @test abs(err) < 1e-14     

    # Right Multiply

    L = PauliOperators.vectorized_rmul(H)

    @test abs(err) < 1e-14     
    a = L*ρv 

    err = norm(Matrix(a.ps) - Matrix(ρmat * Hmat))
    @test abs(err) < 1e-14     

    # Left Multiply

    L = PauliOperators.vectorized_lmul(H)

    @test abs(err) < 1e-14     
    a = L*ρv 

    err = norm(Matrix(a.ps) - Matrix(Hmat * ρmat))
    @test abs(err) < 1e-14     

    display(L)
end
# run()