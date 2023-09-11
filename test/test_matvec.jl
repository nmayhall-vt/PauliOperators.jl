using PauliOperators
using LinearAlgebra
using Random
using Test
using BlockDavidson
# using LinearMaps

@testset "matvec" begin
# function run(N)

    N = 8
    Random.seed!(2)

    a = Pauli(N)

    v0 = rand(ComplexF64,2^N,3)
    v0 = v0 * sqrt(inv(v0'*v0))

    # display(v0'*v0)

    H = random_Pauli(N)
    Hmat = Matrix(H)
    
    @test norm(Hmat*v0 - H*v0) < 1e-14
    
    
    
    a = Pauli(N)
    H = a + a'

    for i in 1:40
        ai = rand()*random_Pauli(N)
        H += ai + ai'
    end
    Hmat = Matrix(H)
    
    @test norm(Hmat*v0 - H*v0) < 1e-14
    @test norm(Hmat*Hmat*v0 - H*H*v0) < 1e-12


    e0, v0 = eigen(Hmat)

    M = 4
    e0 = e0[1:M]
    v0 = v0[:, 1:M]
    for i in 1:M
        display(e0[i])
    end

    dav = Davidson(Hmat, T=ComplexF64, nroots=M)
    @time e1, v1 = eigs(dav)

    @show abs(tr(abs.(v1'*v0)) - M)     
    @test abs(tr(abs.(v1'*v0)) - M) < 1e-4
    
    
    function mymatvec(v) 
        return H*v
    end

    lmat = LinOpMat{ComplexF64}(mymatvec, 2^N, true)
    
    dav = Davidson(lmat, T=ComplexF64, nroots=M)
    @time e2, v2 = eigs(dav)
    
    @show abs(tr(abs.(v2'*v0)) - M)     
    @test abs(tr(abs.(v2'*v0)) - M) < 1e-4
end
# run(4)