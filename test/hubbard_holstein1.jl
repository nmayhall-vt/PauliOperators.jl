using PauliOperators
using LinearAlgebra
using BlockDavidson
import LinearMaps

function hubbard_holstein(;t=1, U=1, ω=.5, g=.2)

   
    N = 4
    Nb = 1      # number of qubits per boson
    pbc = true 

    hopping = [
        (1,2),
        (2,3),
        (3,4),
        (5,6),
        (6,7),
        (7,8)
    ]
    if pbc
        push!(hopping, (4,1))
        push!(hopping, (8,5))
    end

    repulsion = [
        (1,5),
        (2,6),
        (3,7),
        (4,8)
    ]

    Hf = Pauli(2N)*0.0

    for (i,j) in hopping
        
        ai = jordan_wigner(i,2N)'
        aj = jordan_wigner(j,2N)'

        Hf += t * (ai'*aj + aj'*ai)
    end
    for (i,j) in repulsion 
        
        ai = jordan_wigner(i,2N)'
        aj = jordan_wigner(j,2N)'

        Hf += U * (ai'*ai * aj'*aj)
    end

    bi = boson_binary_transformation(Nb)'
    Hb = ω * bi' * bi 
    for i in 2:N
        Hb = Hb ⊕ (ω * bi'*bi)
    end

    H = Hf ⊕ Hb

    aα = jordan_wigner(1, 2N)'
    aβ = jordan_wigner(5, 2N)'
    H += g * (aα' * aα + aβ' * aβ) ⊗ (bi + bi') ⊗ Pauli(3Nb)
    
    aα = jordan_wigner(2, 2N)'
    aβ = jordan_wigner(6, 2N)'
    H += g * (aα' * aα + aβ' * aβ) ⊗ Pauli(1Nb) ⊗ (bi + bi') ⊗ Pauli(2Nb)
    
    aα = jordan_wigner(3, 2N)'
    aβ = jordan_wigner(7, 2N)'
    H += g * (aα' * aα + aβ' * aβ) ⊗ Pauli(2Nb) ⊗ (bi + bi') ⊗ Pauli(1Nb)
    
    aα = jordan_wigner(4, 2N)'
    aβ = jordan_wigner(8, 2N)'
    H += g * (aα' * aα + aβ' * aβ) ⊗ Pauli(3Nb) ⊗ (bi + bi')


    Nop = Pauli(2N)*0.0
    Sz = Pauli(2N)*0.0
    for i in 1:N
        a = jordan_wigner(i, 2N)'
        b = jordan_wigner(i+N, 2N)'
        
        Nop += a'*a + b'*b
        Sz += a'*a + -1 * b'*b
    end
    
    Nop = Nop ⊗ Pauli(N*Nb)
    Sz = Sz ⊗ Pauli(N*Nb)

    clip!(Hb)
    clip!(Hf)
    clip!(H)
    clip!(Nop)
    clip!(Sz)

    return Hf, Hb, H, Nop, Sz
end

Hf, Hb, H, N, Sz = hubbard_holstein()

# e,v = eigen(Matrix(H))
# ni = real.(diag(v'*Matrix(N)*v))
# states = [i ≈ 4 for i in ni];
# display(e[states][1:60])

ndim = 2^nqubits(H)
nvec = 2
function mymatvec(v::Matrix) 
    # display(typeof(v))
    return H * v
end



lmap = LinearMaps.LinearMap(mymatvec, ndim, ndim; issymmetric=false, ismutating=false, ishermitian=true)
lmap = LinOpMat{ComplexF64}(mymatvec, ndim, false)
vin = rand(ComplexF64, ndim, nvec)
# σ = lmap * 
display(typeof(H))
display(typeof(vin))
display(lmap*vin)
# display(lmap * rand(T,ndim, nvec))
dav = Davidson(lmap; max_iter=200, max_ss_vecs=8, tol=1e-6, nroots=nvec, T=ComplexF64)
e, v = eigs(dav)

# dav = Davidson(LinOpMat(H); max_iter=200, max_ss_vecs=8, tol=1e-6, nroots=6, lindep_thresh=1e-10)
# e,v = eigs(dav)
