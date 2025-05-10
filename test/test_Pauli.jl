using PauliOperators
using Test
using LinearAlgebra

@testset "convert" begin
    N = 5
    for i in 1:100
        a = rand(Pauli{N})
        err = Matrix(a) - coeff(a)*Matrix(PauliBasis(a))
        @test norm(err) < 1e-14 
    end
end
