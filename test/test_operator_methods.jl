using PauliOperators
using Test
using Printf
using LinearAlgebra
using Random


@testset "operator_methods" begin
    Random.seed!(1)
  
    N  = 3
    types = []
    push!(types, PauliBasis{N})
    push!(types, Pauli{N})
    push!(types, PauliSum{N, ComplexF64})
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
    types = []
    push!(types, Pauli{N})
    push!(types, PauliSum{N,ComplexF64})
    push!(types, Dyad{N})
    push!(types, DyadSum{N,ComplexF64})

    for T in types
        # Negate
        for i in 1:10
            a = rand(T)
            err = norm(-Matrix(a) - Matrix(-a)) < 1e-9 
            if !err
                @show T, err
            end
            @test err
        end
        
    end
end
        
@testset "Multiplication" begin
    N  = 3
    types = []
    push!(types, PauliBasis{N})
    push!(types, Pauli{N})
    push!(types, PauliSum{N, ComplexF64})
    push!(types, DyadBasis{N})
    push!(types, Dyad{N})
    push!(types, DyadSum{N, ComplexF64})

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
    types1 = []
    push!(types1, PauliBasis{N})
    push!(types1, Pauli{N})
    push!(types1, PauliSum{N, ComplexF64})
    push!(types1, DyadBasis{N})
    push!(types1, Dyad{N})
    push!(types1, DyadSum{N, ComplexF64})
    
    N = 2
    types2 = []
    push!(types2, PauliBasis{N})
    push!(types2, Pauli{N})
    push!(types2, PauliSum{N, ComplexF64})
    push!(types2, DyadBasis{N})
    push!(types2, Dyad{N})
    push!(types2, DyadSum{N, ComplexF64})

    @test Pauli("X") ⊗ Pauli("Y") ⊗ Pauli("Z") == Pauli("XYZ")
    for Ti in 1:length(types1)
        for i in 1:100
            T1 = types1[Ti]
            T2 = types2[Ti]
            a = rand(T1)
            b = rand(T2)
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

@testset "Scalar Mult" begin
    N = 3
    types = []
    push!(types, Pauli{N})
    push!(types, PauliSum{N, ComplexF64})
    push!(types, Dyad{N})
    push!(types, DyadSum{N, ComplexF64})
    for T in types
        for i in 1:100
            # Now scalar multiplication
            
            a = rand(T)
            err = norm(2.3 * Matrix(a)  - Matrix(2.3 * a)) < 1e-14
            if !err
                @show a b err
            end
            @test err
            
            # Now scalar multiplication
            err = norm(2.3 * Matrix(a)'  - Matrix(2.3 * a')) < 1e-14
            if !err
                @show a b err
            end
            @test err
        end
    end
end