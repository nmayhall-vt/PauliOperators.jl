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
   
    # sum
    a = rand(Dyad{N})
    b = rand(Dyad{N})
    err = Matrix(a)+Matrix(b) - Matrix(a+b)
    @test isapprox(norm(err),0, atol=1e-14)
    err = Matrix(a)-Matrix(b) - Matrix(a-b)
    @test isapprox(norm(err),0, atol=1e-14)
    d = rand(Dyad{N}) + rand(Dyad{N})
    e = rand(Dyad{N}) + rand(Dyad{N})
    @show d
    @show e
    println()
    # return
    @show d + e
    @show d - e
    err = Matrix(d)+Matrix(e) - Matrix(d+e)
    @test isapprox(norm(err),0, atol=1e-14)
    err = Matrix(d)-Matrix(e) - Matrix(d-e)
    @test isapprox(norm(err),0, atol=1e-14)
    err = -Matrix(d) - Matrix(-d)
    @test isapprox(norm(err),0, atol=1e-14)
    
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
        err = Matrix(x)*Matrix(s) - Matrix(x*s)
        @test isapprox(norm(err),0, atol=1e-14)
    end
    

    for i in 1:10
        A = PauliSum(N)
        ρ = DyadSum(N)
        for i in 1:4
            A += rand(ScaledPauli{N})
            ρ += rand(ScaledDyad{N,Float64})
        end
        ρ += ρ'
        err = Matrix(A)*Matrix(ρ) - Matrix(A*ρ)
        @test isapprox(norm(err),0, atol=1e-14)
        err = Matrix(A')*Matrix(ρ) - Matrix(A'*ρ)
        @test isapprox(norm(err),0, atol=1e-14)
        err = Matrix(ρ)*Matrix(A) - Matrix(ρ*A)
        @test isapprox(norm(err),0, atol=1e-14)
    end
    for i in 1:10
        A = DyadSum(N)
        B = DyadSum(N)
        for i in 1:4
            A += rand(ScaledDyad{N,Float64})
            B += rand(ScaledDyad{N,Float64})
        end
        err = Matrix(A)*Matrix(B) - Matrix(A*B)
        @test isapprox(norm(err),0, atol=1e-14)
    end
end
# run()