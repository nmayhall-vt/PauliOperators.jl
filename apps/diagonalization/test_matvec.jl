using PauliOperators
using LinearAlgebra
using Random
using Test
# using LinearMaps
using BenchmarkTools
using BlockDavidson
using KrylovKit
using Printf
    
struct LinOp{T} <: AbstractMatrix{T}
    matvec
    dim::Int
    sym::Bool
end

Base.size(lop::LinOp{T}) where {T} = return (lop.dim,lop.dim)
Base.:(*)(lop::LinOp{T}, v::AbstractVector{T}) where {T} = return lop.matvec(v)
Base.:(*)(lop::LinOp{T}, v::AbstractMatrix{T}) where {T} = return lop.matvec(v)
issymmetric(lop::LinOp{T}) where {T} = return lop.sym

function run_diag_test()


    N = 8 
    Random.seed!(2)

    a = Pauli(N)

    v0 = rand(ComplexF64,Int128(2)^N,3)
    v0 = v0 * sqrt(inv(v0'*v0))

    # display(v0'*v0)

    H = rand(Pauli{N})
    Hmat = Matrix(H)
    
    @test norm(Hmat*v0 - H*v0) < 1e-14
    
    
    
    a = Pauli(N)
    H = a + a'

    for i in 1:40
        ai = rand()*rand(Pauli{N})
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
    for i in 1:M
        display(e1[i])
    end

    @show abs(det(v1'*v0))
    @test isapprox(abs(det(v1'*v0)), 1)
    
    scr = deepcopy(v0)
    function mymatvec(v)
        fill!(scr,0.0) 
        mul!(scr, H, v)
        return scr 
    end

    lmat = LinOpMat{ComplexF64}(mymatvec, Int128(2)^N, true)
    
    dav = Davidson(lmat, T=ComplexF64, nroots=M)
    @time e2, v2 = eigs(dav)
    for i in 1:length(e2)
        @printf(" %4i %12.8f %12.8fi\n", i, real(e2[i]), imag(e2[i]))
    end
    
    @show abs(det(v2'*v0))
    @test isapprox(abs(det(v2'*v0)), 1)
    
    println()
    println(" ######## Now test KrylovKit") 
    
    time = @elapsed e3, v3, info = KrylovKit.eigsolve(Hmat, M, :SR, 
    verbosity   = 1, 
    maxiter     = 100, 
    #krylovdim   = max_ss_vecs, 
    issymmetric = true, 
    ishermitian = true, 
    eager       = true,
    tol         = 1e-5)
    
    for i in 1:length(e3)
        @printf(" %4i %12.8f %12.8fi\n", i, real(e3[i]), imag(e3[i]))
    end


    # The following is not yet working!


    # lmat = LinOp{ComplexF64}(mymatvec, Int128(2)^N, true)

    # randmat = rand(size(lmat)...)
    # print(size(Hmat))
    # print(size(randmat))
    # print(size(lmat))
    # println(typeof(Hmat))
    # println(typeof(lmat*randmat))
    # @test isapprox(norm(Matrix(lmat*randmat) - Hmat*randmat), 0)

    # time = @elapsed e3, v3, info = KrylovKit.eigsolve(lmat, M, :SR, 
    # verbosity   = 1, 
    # maxiter     = 100, 
    # #krylovdim   = max_ss_vecs, 
    # issymmetric = true, 
    # ishermitian = true, 
    # eager       = true,
    # tol         = 1e-5)
    
    # for i in 1:length(e3)
    #     @printf(" %4i %12.8f %12.8fi\n", i, real(e3[i]), imag(e3[i]))
    # end

end
run_diag_test()