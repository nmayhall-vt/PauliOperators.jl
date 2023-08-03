using PauliOperators 
using Test
# using BenchmarkTools 
# using Printf
# using LinearAlgebra
# using Random

@testset "PauliBitString" begin

    println() 
    a = PauliBitString("XYZIXY")
    b = PauliBitString("YXXYZZ")
    c = PauliBitString("ZZYYYX")
    c = rotate_phase(c,1)
    display(a)
    display(b)
    display(a*b)
    display(c)
    @test c == a*b 

    @test commute(a,b) == false
    @test get_phase(c) == -1
    @test get_phase(negate(c)) == 1
    @test get_phase(rotate_phase(c,0)) == -1
    @test get_phase(rotate_phase(c,1)) == -1im
    @test get_phase(rotate_phase(c,2)) == 1
    @test get_phase(rotate_phase(c,3)) == 1im
    @test get_phase(rotate_phase(c,5)) == -1im



    ZX = [0+0im  1+0im   0+0im   0+0im
    1+0im  0+0im   0+0im   0+0im
    0+0im  0+0im   0+0im  -1+0im
    0+0im  0+0im  -1+0im   0+0im]

    @test all(ZX .== Matrix(PauliBitString("ZX")))
end