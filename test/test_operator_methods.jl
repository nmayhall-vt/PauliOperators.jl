using PauliOperators
using Test
using Printf
using LinearAlgebra
using Random


@testset "operator_methods" begin
    Random.seed!(1)
  
    N  = 3
    types = [Pauli{N}, Dyad{N}, PauliSum{N, ComplexF64}]
    push!(types, Dyad{N})
    push!(types, DyadBasis{N})
    push!(types, DyadSum{N, ComplexF64})

    for T in types
        # Hermiticity
        for i in 1:10
            a = rand(T)
            if ishermitian(a)
                @test norm(Matrix(a) - Matrix(a)') < 1e-14
            end
        end 
    
        # Adjoint
        for i in 1:10
            a = rand(T)
            err = norm(Matrix(a)' - Matrix(a')) < 1e-9 
            if !err
                @show a err
            end
            @test err
        end
        
        # Trace 
        for i in 1:10
            a = rand(T)
            err = abs(tr(Matrix(a)) - tr(a)) < 1e-9 
            if !err
                @show a, err
            end
            @test err
        end
        
        # Addition
        # for a in T 
        #     for b in T 
        #         err = norm(Matrix(a) + Matrix(b) - Matrix(a+b)) < 1e-9
        #         if !err
        #             @show a, b, err
        #         end
        #         @test err
        #     end
        # end
        for i in 1:100
            a = rand(T)
            b = rand(T)
            err = norm(Matrix(a) + Matrix(b) - Matrix(a+b)) < 1e-9
            if !err
                @show T, err
            end
            @test err
            
            err = norm(Matrix(a) + Matrix(b)' - Matrix(a+b')) < 1e-9
            if !err
                @show T, err
            end
            @test err
            
            err = norm(Matrix(a)' + Matrix(b) - Matrix(a'+b)) < 1e-9
            if !err
                @show T, err
            end
            @test err
        end
        
    end
end
@testset "Negate" begin
    Random.seed!(1)
  
    N  = 3
    types = [Pauli{N}, PauliSum{N, ComplexF64}]

    for T in types
        # Negate
        for i in 1:10
            a = rand(T)
            @test norm(-Matrix(a) - Matrix(-a)) < 1e-9 
        end
        
    end
end
        
@testset "Multiplication" begin
    N  = 3
    types = [Pauli{N}, PauliBasis{N}, PauliSum{N, ComplexF64}, Dyad{N}, DyadBasis{N}, DyadSum{N, ComplexF64}]

    for T in types
        for i in 1:20
            a = rand(T)
            b = rand(T)
            err = norm(Matrix(a) * Matrix(b) - Matrix(a*b)) < 1e-14
            if !err
                @show a b err
            end
            @test err
            
            # Now adjoint multiplication
            err = norm(Matrix(a)' * Matrix(b) - Matrix(a'*b)) < 1e-14
            if !err
                @show a b err
            end
            @test err
        end 
        
    end
end  

@testset "Tensor product" begin
    N = 3
    types = [Pauli{N}, PauliBasis{N}, PauliSum{N, ComplexF64}]
    # push!(types, Dyad{N})
    # push!(types, DyadBasis{N})

    @test Pauli("X") ⊗ Pauli("Y") ⊗ Pauli("Z") == Pauli("XYZ")
    for T in types
        for i in 1:100
            a = rand(T)
            b = rand(T)
            err = norm(kron(Matrix(b), Matrix(a)) - Matrix(a⊗b)) < 1e-14
            if !err
                @show T, err
            end
            @test err
        end 
        
    end
end  

@testset "PauliSum" begin
    Random.seed!(1)
    N = 5
    # Hermiticity
    for i in 1:20
        a = rand(Pauli{N})
        b = rand(Pauli{N})
        @test norm(Matrix(a) + Matrix(b) - Matrix(a+b)) < 1e-14
    end 
    
end  
