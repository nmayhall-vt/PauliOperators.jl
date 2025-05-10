using PauliOperators
using Test
using LinearAlgebra

@testset "Ket" begin
    Random.seed!(1)
    N = 5
    for i in 1:100
        a = rand(Ket{7})
        b = rand(Ket{7})
        err = Vector(a) + Vector(b) - Vector(a+b) 
        @test norm(err) < 1e-14 
       

        err = Vector(a') + Vector(b') - Vector(a'+b') 
        @test norm(err) < 1e-14 
    end
end
