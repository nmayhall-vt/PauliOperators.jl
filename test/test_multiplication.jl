using PauliOperators
using Test
using Printf
using LinearAlgebra
using Random

@testset "Multiplication" begin
    Random.seed!(1)
  
    N  = 3
    types = []
    push!(types, PauliBasis{N})
    push!(types, Pauli{N})
    # push!(types, PauliSum{N, ComplexF64})
    push!(types, DyadBasis{N})
    push!(types, Dyad{N})
    # push!(types, DyadSum{N, ComplexF64})

    for T1 in types
        for T2 in types
            for i in 1:10
                a = rand(T1)
                b = rand(T2)
                err = norm(Matrix(a)*Matrix(b) - Matrix(a*b)) < 1e-14
                if err == false
                    @show a b err
                end
                @test err
            end 
        end 
    end 
end 
    
