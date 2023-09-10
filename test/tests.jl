using PauliOperators 
using Test
# using BenchmarkTools 
# using Printf
using LinearAlgebra
using Random

@testset "Pauli" begin
    Random.seed!(1)

    println() 
    a = Pauli("XYZIXY")
    b = Pauli("YXXYZZ")
    c = Pauli("ZZYYYX")
    c = rotate_phase(c,1)
    # display(a)
    # display(b)
    # display(a*b)
    # display(c)
    @test c == a*b 
  
    println(" Now PF")
    a = FixedPhasePauli("XYZIXY")
    b = FixedPhasePauli("YXXYZZ")
    display(a)
    display(b)
    display(a*b)
    display(get_phase(a,b)*a*b)
    display(1*c)
    display(c)
    @test c == Pauli{6}(PauliOperators.phase(a,b), a*b) 
    @test 1*c == get_phase(a,b)*a*b 
    @test 2*c == 2*a*b*get_phase(a,b) 

    @test commute(a,b) == false
    @test get_phase(c) == -1
    @test get_phase(negate(c)) == 1
    @test get_phase(rotate_phase(c,0)) == -1
    @test get_phase(rotate_phase(c,1)) == -1im
    @test get_phase(rotate_phase(c,2)) == 1
    @test get_phase(rotate_phase(c,3)) == 1im
    @test get_phase(rotate_phase(c,5)) == -1im

    a = FixedPhasePauli("ZXYI"); b = FixedPhasePauli("YZXX");
    @test norm(Matrix(a*b)*get_phase(a,b) - Matrix(a)*Matrix(b)) ≈ 0

    display((a.z, a.x))
    display((b.z, b.x))
    display((c.p.z, c.p.x))

    println()
    a = Pauli("XYZIXY")
    b = Pauli("YXXYZZ")
    c = Pauli("ZZYYYX")
    c = rotate_phase(c,1)
    @test a < b
    @test b > a
    @test a > c
    @test c < a
    @test b > c
    @test c < b
    @test c == c
  
    # make sure that two operators that only differ in phase are equal in hash
    d = rotate_phase(a,1)
    # @test hash(phasefree(a)) == hash(phasefree(d))
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
    for i in 1:10
        N = 8
        a = random_Pauli(N)
        b = random_Pauli(N)
        c = random_Pauli(N)
        d = random_Pauli(N)

        s1 = a + b + d
        s2 = a + c + d
        display(s1 * s2)
        # println()
        # display(s2)

        @test all(Matrix(s1) * Matrix(s2) - Matrix(s1 * s2) .≈ 0)
    end


    a = random_Pauli(6)
    b = random_Pauli(6)
    s = a + b
    c = 1.23
    @test all(c*Matrix(s) - Matrix(c * s) .≈ 0)
    println()

    ZX = [0+0im  1+0im   0+0im   0+0im
    1+0im  0+0im   0+0im   0+0im
    0+0im  0+0im   0+0im  -1+0im
    0+0im  0+0im  -1+0im   0+0im]

    @test all(ZX .== Matrix(Pauli("ZX")))


    bdag = boson_binary_transformation(6)
    evals = real.(eigvals(Matrix(bdag*adjoint(bdag))))
    @test all(abs.(evals .- range(0,63)) .< 1e-12)

    println("Print JW Transformation")
    fd = PauliOperators.jordan_wigner(5,10)
    tmp = .5*(Pauli("ZZZZXIIIII") + -1im*Pauli("ZZZZYIIIII"))
    @test all(Matrix(fd) .≈ Matrix(tmp))
    @test fd ≈ tmp


    # Tensor products
    @test Pauli("X") ⊗ Pauli("Y") ⊗ Pauli("Z") == Pauli("XYZ")
    
    @test Pauli("X") ⊗ (Pauli("Y") + Pauli("Z")) ≈ Pauli("XY") + Pauli("XZ") 
    
    @test (Pauli("X") + Pauli("Y")) ⊗ (Pauli("Y") + Pauli("Z")) ≈ Pauli("XY") + Pauli("XZ") + Pauli("YY") + Pauli("YZ") 


    # is_hermitian
    @test is_hermitian(Pauli("XXXXX")) == true
    @test is_hermitian(Pauli("XXXXY")) == true
    @test is_hermitian(Pauli("XXYZY")) == true
    @test is_hermitian(Pauli("XXXXX")) == true
    @test is_hermitian(phasefree(Pauli("YXXXX"))) == false 
    @test is_hermitian(phasefree(Pauli("YYYXX"))) == false 
    @test is_hermitian(phasefree(Pauli("YZYXX"))) == true 

    # Direct sums
    @test Pauli("X") ⊕ Pauli("IZ") ≈ Pauli("XII") + Pauli("IIZ")
    @test Pauli("IZ") ⊕ (Pauli("Y") + Pauli("Z")) ≈ Pauli("IZI") + Pauli("IIY") + Pauli("IIZ")
    @test (Pauli("YI") + Pauli("ZX")) ⊕ Pauli("IZ") ≈ Pauli("YIII") + Pauli("ZXII") + Pauli("IIIZ")
    @test (Pauli("YI") + Pauli("ZX")) ⊕ (Pauli("Z") + Pauli("X")) ≈ Pauli("YII") + Pauli("ZXI") + Pauli("IIZ") + Pauli("IIX")
    
    # Pauli KetBitString multiplication 
    @test Pauli("YZYZ") * KetBitString([1,0,0,0]) == (1, KetBitString([0,0,1,0]))
    @test Pauli("YZYZ") * KetBitString([1,0,1,0]) == (-1, KetBitString([0,0,0,0]))
    @test Pauli("YZYZ") * KetBitString([1,1,1,0]) == (1, KetBitString([0,1,0,0]))
    @test Pauli("YZXZ") * KetBitString([1,1,1,0]) == (1im, KetBitString([0,1,0,0]))
    @test Pauli("XZXZ") * KetBitString([1,1,1,0]) == (-1, KetBitString([0,1,0,0]))

    # ScaledPauli
    N=8
    for i in 1:10
        a = ScaledPauli(random_Pauli(N))
        b = ScaledPauli(random_Pauli(N))

        a *= 2.3
        b *= 3.2

        @test get_coeff(a * b) ≈ get_coeff(a) * get_coeff(b) * get_phase(a.pauli, b.pauli)
    end

    # Test unique!
    # Create vector of scaled paulis
    v = [ScaledPauli(random_Pauli(8)) for i in 1:10]
   
    # add duplicates
    w = deepcopy(v)
    for i in 1:10
        for j in i:10
            push!(w, v[rand(1:10)]*i)
        end
    end
    display(v)
    println() 
    display(w)

    @test length(v) != length(w)
    @test length(v) == length(unique(w))
    
    mat1 = Matrix(w)
    mat2 = Matrix(unique(w))

    @test norm(mat1-mat2) ≈ 0
end
