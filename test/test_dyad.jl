using PauliOperators
using LinearAlgebra
using Random
using Test
using BenchmarkTools

@testset "dyad" begin
# function run()
    N = 8 
    Random.seed!(2)

    a = Pauli(N)
    di = Dyad([0],[0])
    display(di)
    di = rand(Dyad{N}) 
    dj = rand(Dyad{N}) 
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
    @test dl.dyad.ket == di.ket 
    @test dl.dyad.bra == dk.dyad.bra

    @show di'*dj 
    @show di*dj 
    @test (di*di').dyad == Dyad(di.ket, adjoint(di.ket))


    # Now test Pauli*Dyad multiplication
    x = Pauli("X")
    s = ScaledDyad(1, 1, 0, 0)
    @show x
    @show s
    @show x.pauli * s.dyad
    @show x * s.dyad
    @show (2.0*x) * s
    # @show typeof((x * s.dyad.ket)[1])

    println()
    x = Pauli("Y")
    s = ScaledDyad(1, 1, 0, 0)
    @show x.pauli
    @show s.dyad
    @show x.pauli * s.dyad
    @show s.dyad * x.pauli
   
    # return
    N = 8
    for i in 1:10
        x = rand(ScaledPauli{N})
        s = rand(ScaledDyad{N})
        err = Matrix(x)*Matrix(s) - Matrix(x*s)
        @test isapprox(norm(err),0, atol=1e-14)
    end
    
    for i in 1:10
        x = rand(Pauli{N})
        s = rand(ScaledDyad{N})
        err = Matrix(s)*Matrix(x) - Matrix(s*x)
        @test isapprox(norm(err),0, atol=1e-14)
    end
    
    for i in 1:10
        x = rand(ScaledPauli{N})
        s = rand(ScaledDyad{N})
        err = Matrix(s)*Matrix(x) - Matrix(s*x)
        @test isapprox(norm(err),0, atol=1e-14)
    end
    
    # for i in 1:10
    #     p = rand(ScaledPauli{N})
    #     k = rand(KetBitString{N})
    #     b = rand(BraBitString{N})
    #     err = Matrix(b)*Matrix(p)*Matrix(k) - Matrix(b*p*k)
    #     @test isapprox(norm(err),0, atol=1e-14)
    # end
end
# run()