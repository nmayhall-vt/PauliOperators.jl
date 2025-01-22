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
    @test c == a*b 
    @test norm(-Matrix(a) - Matrix(-a)) < 1e-9 
  
    println(" Now PF")
    a = FixedPhasePauli("XYZIXY")
    b = FixedPhasePauli("YXXYZZ")

    @test commute(a,b) == false
    @test get_phase(c) == -1
    @test get_phase(negate(c)) == 1
    @test get_phase(rotate_phase(c,0)) == -1
    @test get_phase(rotate_phase(c,1)) == -1im
    @test get_phase(rotate_phase(c,2)) == 1
    @test get_phase(rotate_phase(c,3)) == 1im
    @test get_phase(rotate_phase(c,5)) == -1im

    a = FixedPhasePauli("ZXYI"); 
    b = FixedPhasePauli("YZXX");
    @test norm(Matrix(a*b) * get_phase(a,b) - Matrix(a)*Matrix(b)) ≈ 0

    display((a.z, a.x))
    display((b.z, b.x))
    display((c.pauli.z, c.pauli.x))

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
   

    ZX = [0+0im  0+0im   1+0im   0+0im
          0+0im  0+0im   0+0im  -1+0im
          1+0im  0+0im   0+0im   0+0im
          0+0im -1+0im   0+0im   0+0im]

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


    # Test unique!
    # Create vector of scaled paulis
    v = [ScaledPauli(rand(Pauli{8})) for i in 1:10]
   
    # add duplicates
    w = deepcopy(v)
    for i in 1:10
        for j in i:10
            push!(w, v[rand(1:10)]*i)
        end
    end
    # display(v)
    # println() 
    # display(w)

    @test length(v) != length(w)
    @test length(v) == length(unique(w))
    
    mat1 = Matrix(w)
    mat2 = Matrix(unique(w))

    @test norm(mat1-mat2) ≈ 0

    ##      test matvec
    N = 8
    H = rand(Pauli{N})
    for i in 1:40
        H += rand()*rand(Pauli{N})
    end
   
    # Vector as Matrix
    v = rand(ComplexF64, Int128(2)^N,1)
    v .= v ./ norm(v) 
    @test norm(Matrix(H)*v - H*v) < 1e-13 
    
    # Matrix
    v = rand(ComplexF64, Int128(2)^N,3)
    v = v *  sqrt(inv(v'*v)) 
    @test norm(Matrix(H)*v - H*v) < 1e-13 
    
    
    # Vector 
    v = rand(ComplexF64, Int128(2)^N)
    v .= v ./ norm(v) 
    @show norm(Matrix(H)*v - H*v) 
    @test norm(Matrix(H)*v - H*v) < 1e-13 

    
    ##      test matvec with a scaled pauli vector
    N = 8
    H = ScaledPauliVector(N) 
    for i in 1:40
        push!(H, rand(ScaledPauli{N}))
    end
    @show norm(Matrix(H)*v - H*v) 
    @test norm(Matrix(H)*v - H*v) < 1e-13 


    ##      diag
    H = PauliSum(N)
    for i in 1:100
        sum!(H, rand(ScaledPauli{N}))
    end
    @test abs(tr(H) - tr(Matrix(H))) < 1e-9

    test = true
    for (op,coeff) in diag(H).ops
        test = test && is_diagonal(op)
    end
    H = PauliSum(N)
    for i in 1:100
        sum!(H, FixedPhasePauli{N}(rand(0:Int128(2)^N-1),0))
    end
    @test length(diag(H)) == length(H)


    ##      test SparseKetBasis
    for i in 1:10
        v = SparseKetBasis(8)
        for i in 1:2
            sum!(v, rand(KetBitString{8}), rand())
        end
        # @show dot(v,v) - Vector(v)'*Vector(v)  
        @test dot(v,v) - Vector(v)'*Vector(v) < 1e-14 
  
        a = rand()
        @test norm(a*Vector(v) - Vector(a*v)) < 1e-4
        
        o = rand(Pauli{8})
        # @show norm(Matrix(o)*Vector(v) - Vector(o*v)) 
        @test norm(Matrix(o)*Vector(v) - Vector(o*v)) < 1e-13 

        o = rand(FixedPhasePauli{8})
        # @show norm(Matrix(o)*Vector(v) - Vector(o*v)) 
        @test norm(Matrix(o)*Vector(v) - Vector(o*v)) < 1e-13 

        o = rand(ScaledPauli{8})
        # @show norm(Matrix(o)*Vector(v) - Vector(o*v)) 
        @test norm(Matrix(o)*Vector(v) - Vector(o*v)) < 1e-13 

        o = rand(ScaledPauli{8}) + rand(ScaledPauli{8})
        for i in 1:10
            o += rand(ScaledPauli{8})
        end
        # @show norm(Matrix(o)*Vector(v) - Vector(o*v)) 
        @test norm(Matrix(o)*Vector(v) - Vector(o*v)) < 1e-13 
    end
    
    
    ##      sum/Subtract 
    N = 8
    H1 = PauliSum(N)
    H2 = PauliSum(N)
    for i in 1:100
        sum!(H1, rand(ScaledPauli{N}))
        sum!(H2, rand(ScaledPauli{N}))
    end
    @test norm(Matrix(H1) + Matrix(H2) - Matrix(H1 + H2)) < 1e-8
    @test norm(Matrix(H1) - Matrix(H2) - Matrix(H1 - H2)) < 1e-8

    ket = rand(KetBitString{N})
    @test abs(expectation_value(H1, ket) - Vector(ket)' * Matrix(H1)*Vector(ket)) < 1e-8
    @test abs(expectation_value(H2, ket) - Vector(ket)' * Matrix(H2)*Vector(ket)) < 1e-8

    # KetBitString indexing
    @test string(KetBitString(4,0)) == "|0000>"
    @test string(KetBitString(4,1)) == "|1000>"
    @test string(KetBitString(4,2)) == "|0100>"
    @test string(KetBitString(4,3)) == "|1100>"
    @test string(KetBitString(4,4)) == "|0010>"
    @test string(KetBitString(4,5)) == "|1010>"
    @test string(KetBitString(4,6)) == "|0110>"
    @test string(KetBitString(4,7)) == "|1110>"
    @test string(KetBitString(4,8)) == "|0001>"
    @test string(KetBitString(4,9)) == "|1001>"
    @test string(KetBitString(4,10)) == "|0101>"
    @test string(KetBitString(4,11)) == "|1101>"
    @test string(KetBitString(4,12)) == "|0011>"
    @test string(KetBitString(4,13)) == "|1011>"
    @test string(KetBitString(4,14)) == "|0111>"
    @test string(KetBitString(4,15)) == "|1111>"

    println("Test Commutator")
    spv1 = (1/2) .* [
                    ScaledPauli(Pauli(4; X=[1], Y=[3])),
                   -ScaledPauli(Pauli(4; X=[3], Y=[1]))];
    spv2 = (1/2) .* [
                    ScaledPauli(Pauli(4; X=[2], Y=[4])),
                   -ScaledPauli(Pauli(4; X=[4], Y=[2]))];
   
    display(spv1)
    println()
    display(spv2)
    println()
    spv_comm = commutator(spv1,spv2)
    display(spv1)
    println()
    display(spv2)
    println()
    display(spv_comm)
    return
    @test isempty(spv_comm)

    sp1 = ScaledPauli(Pauli("X")); sp2 = ScaledPauli(Pauli("Y"));
    sp_comm = commutator(sp1,sp2)
    sp_comm_true = (0.0+2.0im)*[ScaledPauli(Pauli("Z"))]
    @test sp_comm == sp_comm_true



    println("Test Adjoint")
    A = PauliSum(5)
    B = PauliSum(5)
    for i in 1:50
        sum!(A, rand(ScaledPauli{5}))
        sum!(B, rand(ScaledPauli{5}))
    end
    @test norm(Matrix(A)' - Matrix(A')) < 1e-8
    @test norm(Matrix(B)' - Matrix(B')) < 1e-8
end




@testset "Mul1" begin
    # Test Multiply
    ntests = 10 
    N = 7
    types = [Pauli{N}, FixedPhasePauli{N}, ScaledPauli{N}, PauliSum{N}, Dyad{N}, ScaledDyad{N}, DyadSum{N}]
    for t1 in types
        for t2 in types
            for i in 1:ntests
                a = rand(t1)
                b = rand(t2)
                phase = 1
                if t1 == FixedPhasePauli{N} && t2 == FixedPhasePauli{N}
                    phase = get_phase(a,b)
                end
                @test all(abs.(Matrix(a) * Matrix(b) - Matrix(a * b) * phase) .< 1e-14)
            end
        end
    end
end
@testset "Add1" begin
    # Test Multiply
    ntests = 10 
    N = 7
    types = [Pauli{N}, FixedPhasePauli{N}, ScaledPauli{N}, PauliSum{N}]
    for t1 in types
        for t2 in types
            for i in 1:ntests
                a = rand(t1)
                b = rand(t2)
                @test all(abs.(Matrix(a) + Matrix(b) - Matrix(a + b)) .< 1e-14)
            end
        end
    end
end
@testset "Add2" begin
    # Test Multiply
    ntests = 10 
    N = 7
    types = [Dyad{N}, ScaledDyad{N}]
    for t1 in types
        for t2 in types
            for i in 1:ntests
                a = rand(t1)
                b = rand(t2)
                @test all(abs.(Matrix(a) + Matrix(b) - Matrix(a + b)) .< 1e-14)
            end
        end
    end
end

@testset "is_hermitian" begin
    # is_hermitian
    @test is_hermitian(Pauli("XXXXX")) == true
    @test is_hermitian(Pauli("XXXXY")) == true
    @test is_hermitian(Pauli("XXYZY")) == true
    @test is_hermitian(Pauli("XXXXX")) == true
    @test is_hermitian(phasefree(Pauli("YXXXX"))) == false 
    @test is_hermitian(phasefree(Pauli("YYYXX"))) == false 
    @test is_hermitian(phasefree(Pauli("YZYXX"))) == true 
end
@testset "expectation_value" begin
    N = 8
    types = [Pauli{N}, FixedPhasePauli{N}, ScaledPauli{N}, PauliSum{N}]
    for ti in types
        for i in 1:10
            o = rand(ti)
            v = rand(KetBitString{N})
            @test expectation_value(o, v) == Vector(v)'*Matrix(o)*Vector(v) 
        end
    end
end