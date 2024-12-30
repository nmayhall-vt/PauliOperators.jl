using PauliOperators
using LinearAlgebra
using Printf

function build_full_hamiltonian(ω,λ,N)
    "https://www.uni-ulm.de/fileadmin/website_uni_ulm/nawi.inst.260/%C3%9Cbungen/6Spin_Boson_Modell.pdf"
    Hs = PauliSum(1)
    Hb = PauliSum(N)
    bi = boson_binary_transformation(N)'

    Hs += ω/2 * Pauli(1,Z=[1]) 
    Hb = ω * bi' * bi 
    H = Hs ⊕ Hb

    Z = PauliSum(Pauli("Z"))

    # σ⁺ = .5*(Pauli("X") + 1im*Pauli("Y"))
    # Hc = σ⁺⊗bi
    # Hc += Hc'
    Hc = λ * Z ⊗ (bi'+bi)
    clip!(Hc)
    # println("Hc")
    # display(Hc)
    H += Hc 
    return H
end

function run()
    ω0 = 1
    N = 4
    g = .1 
    H = build_full_hamiltonian(ω0, g, N)

    # println("H")
    # display(H)

    e,v = eigen(Matrix(H))
    nstates = min(length(e),5)
    println(" Eigenvalues")
    for i in 1:nstates
        @printf(" %3i %12.8f + %12.8fi\n",i,real(e[i]), imag(e[i]))
    end

    for i in 1:nstates
        println("\n State ", i)
        v0 = v[:,i]
        v0 = reshape(v0,(2,2^N))
        ρ = v0*v0'
        Z = Matrix(Pauli("Z"))
        println(" ρ:")
        ρ = round.(real(ρ),digits=3)
        
        display(round.(real(ρ),digits=3))

        sz = tr(ρ*Z)
        @printf(" tr(Sz,ρ) = %12.8f %12.8fi\n", real(sz), imag(sz))
    end
end
run()