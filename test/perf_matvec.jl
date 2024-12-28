using PauliOperators
using LinearAlgebra
using Random
using Test
using BlockDavidson
# using LinearMaps
using BenchmarkTools

function run(nqubits, nvecs, nops)

    N = nqubits

    Random.seed!(2)

    v0 = rand(ComplexF64,Int128(2)^N,nvecs)
    v0 = v0 * sqrt(inv(v0'*v0))

    a = Pauli(N)
    H = a + a'

    for i in 1:nops
        ai = rand()*random_Pauli(N)
        H += ai + ai'
    end
    
    # scr = deepcopy(v0)
    # fill!(scr,0.0) 
    # mul!(scr, H, v0)
    # @test norm(Matrix(H)*v0 - scr) < 1e-12
    
    # fill!(scr,0.0) 
    # PauliOperators.mul2!(scr, H, v0)
    # @test norm(Matrix(H)*v0 - scr) < 1e-12
    
    scr = deepcopy(v0)
    
    fill!(scr,0.0) 
    @btime mul!($scr, $H, $v0)
    
    # fill!(scr,0.0) 
    # @btime PauliOperators.mul2!($scr, $H, $v0)
    
   
    v0 = v0[:,1]
    scr = deepcopy(v0)
    fill!(scr,0.0) 
    mul!(scr, H, v0)
    # @test norm(Matrix(H)*v0 - scr) < 1e-12
    @btime mul!($scr, $H, $v0)
end


run(8, 2, 20)