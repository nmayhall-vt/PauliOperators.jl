using PauliOperators
using LinearAlgebra
using Random
using Test
using BenchmarkTools

# @testset "dyad" begin
function run()
    N = 8 
    Random.seed!(2)

    a = Pauli(N)
    di = Dyad([0],[0])
    display(di)
    di = random_Dyad(8) 
    dj = random_Dyad(8) 
    @show di
    @show dj

    @show di+dj
    
    @show di==dj
    @show isequal(di,dj)
    
    @show di-dj
    
    @show di*2
    @show di*2.2
    dk = dj*di
    @show di
    @show dj
    @show dk
    dl = di*dj*dk
    println(Int(dl.dyad.ket.v))
    println(Int(dl.dyad.bra.v))
    @test dl.dyad.ket.v == 9
    @test dl.dyad.bra.v == 57 

    @show di'*dj 
    @show di*dj 
    @show di*di'


    # Now test Pauli*Dyad multiplication
    x = Pauli("X")
    s = ScaledDyad(1, 1, 0, 0)
    @show x
    @show s
    @show x.pauli * s.dyad
    @show x * s.dyad
    @show (2.0*x) * s
    # @show typeof((x * s.dyad.ket)[1])

    x = Matrix(x)
    s = Matrix(s)
    
    N = 8
    for i in 1:10
        x = random_ScaledPauli(N)
        s = random_ScaledDyad(N)
        err = Matrix(x)*Matrix(s) - Matrix(x*s)
        @test isapprox(norm(err),0)
    end
    
    for i in 1:10
        x = random_Pauli(N)
        s = random_ScaledDyad(N)
        err = Matrix(x)*Matrix(s) - Matrix(x*s)
        @test isapprox(norm(err),0)
    end
end
run()