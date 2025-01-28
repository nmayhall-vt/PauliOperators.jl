using PauliOperators
using Test
using LinearAlgebra

@testset "convert" begin
    N = 5
    for i in 1:100
        a = rand(Pauli{N})
        err = Matrix(a) - 1im^global_phase(a)*Matrix(PauliBasis(a))
        @test norm(err) < 1e-14 
       
        b = Pauli{N}(symplectic_phase(a), a.z, a.x)

        err = Matrix(b) - Matrix(Pauli(PauliBasis(a)))
        @test norm(err) < 1e-14 
    end
end
