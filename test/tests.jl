using PauliOperators 
using Test
# using BenchmarkTools 
# using Printf
# using LinearAlgebra
# using Random

@testset "Pauli" begin

    println() 
    a = Pauli("XYZIXY")
    b = Pauli("YXXYZZ")
    c = Pauli("ZZYYYX")
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


    display((a.z, a.x))
    display((b.z, b.x))
    display((c.z, c.x))

    println()
    @test a < b
    @test b > a
    @test a > c
    @test c < a
    @test b > c
    @test c < b
    @test c == c
  
    # make sure that two operators that only differ in phase are equal in hash
    d = rotate_phase(a,1)
    @test hash(a) == hash(d)
    @test a != d
   
    sum1 = a + b
    display(sum1)
    println()

    sum2 = a + d
    display(sum2)

    println()
    sum3 = sum1 + sum2
    sum!(sum1, sum2)
    display(sum1)
    
    println()
    display(sum3)

    check = true
    for key in keys(sum3)
        check = sum3[key] ≈ sum1[key] && check
    end
    @test check

    # Test Multiply
    println("Test Multiply")
    a = Pauli(123, 2345, 5)
    b = Pauli(533, 345, 5)
    c = Pauli(33, 435, 5)
    d = Pauli(513, 534, 5)

    s1 = a + b
    s2 = a + c + d 
    display(s1*s2)
    println()
    display(s2)



    println()
    ZX = [0+0im  1+0im   0+0im   0+0im
    1+0im  0+0im   0+0im   0+0im
    0+0im  0+0im   0+0im  -1+0im
    0+0im  0+0im  -1+0im   0+0im]

    @test all(ZX .== Matrix(Pauli("ZX")))


end