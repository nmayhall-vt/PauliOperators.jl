using PauliOperators
using LinearAlgebra
using Random
using Test
# using LinearMaps
using BenchmarkTools

@testset "matvec" begin
# function run()

    N = 8 
    Random.seed!(2)

    a = Pauli(N)

    v0 = rand(ComplexF64,Int128(2)^N,3)
    v0 = v0 * sqrt(inv(v0'*v0))

    # display(v0'*v0)

    H = rand(Pauli{N})
    Hmat = Matrix(H)
    
    @test norm(Hmat*v0 - H*v0) < 1e-14
    
    
    
    H = PauliSum(N) 
    for i in 1:40
        ai = rand(ScaledPauli{N})
        H += ai + ai'
    end
    Hmat = Matrix(H)
    
    @test norm(Hmat*v0 - H*v0) < 1e-14
    @test norm(Hmat*Hmat*v0 - H*H*v0) < 1e-12
    
    

    ####
    ## Now confirm matvec works as a function
    vrand = rand(ComplexF64,size(Hmat)[1], 4)
    
    function mymatvec(v)
        scr = deepcopy(v)
        fill!(scr,0.0) 
        mul!(scr, H, v)
        return scr 
    end
    
    @test norm(Hmat*vrand - Matrix(mymatvec(vrand))) < 1e-13

   
end
# run()