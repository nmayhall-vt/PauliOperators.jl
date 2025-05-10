using PauliOperators
using Test
using Printf
using LinearAlgebra
using Random

@testset "Allocations" begin
# function test()
    Random.seed!(1)
  
    N  = 9

    i = @allocated for i in 1:20
        a = rand(PauliBasis{N})
        b = rand(Pauli{N})
        c = rand(DyadBasis{N})
        d = rand(Dyad{N})
    end
    @test i==0
end 

