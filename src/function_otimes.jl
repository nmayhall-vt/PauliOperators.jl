otimes(p1::Pauli{N}, p2::Pauli{M}) where {N,M} = Pauli{N+M}((p1.θ + p2.θ)%4, p1.pauli ⊗ p2.pauli)
otimes(p1::FixedPhasePauli{N}, p2::FixedPhasePauli{M}) where {N,M} = FixedPhasePauli{N+M}(p1.z | p2.z << N, p1.x | p2.x << N)
otimes(p::Pauli{N}, fpp::FixedPhasePauli{M}) where {N,M} = otimes(p, Pauli{M}(0,fpp))
otimes(fpp::FixedPhasePauli{M}, p::Pauli{N}) where {N,M} = otimes(Pauli{M}(0,fpp), p)

function otimes(sp::ScaledPauli{N}, p::Pauli{M}) where {N,M}
    out_pauli = sp.pauli ⊗ p.pauli
    out_coeff = sp.coeff * get_phase(p)
    out = ScaledPauli{N+M}(out_coeff, out_pauli)
    return out 
end

function otimes(p::Pauli{M}, sp::ScaledPauli{N}) where {M,N}
    out_pauli = p.pauli ⊗ sp.pauli
    out_coeff = sp.coeff * get_phase(p)
    out = ScaledPauli{N+M}(out_coeff, out_pauli)
    return out 
end

function otimes(ps::PauliSum{N}, p::Pauli{M}) where {N,M}
    out = PauliSum(N+M)
    for (op,coeff) in ps.ops
        out.ops[op ⊗ p.pauli] = coeff * get_phase(p)
    end
    return out 
end
function otimes(p::Pauli{N}, ps::PauliSum{M}) where {N,M}
    out = PauliSum(N+M)
    for (op,coeff) in ps.ops
        out.ops[p.pauli ⊗ op] = coeff * get_phase(p)
    end
    return out 
end

function otimes(p1::PauliSum{N}, p2::PauliSum{M}) where {N,M}
    out = PauliSum(N+M)
    for (op1,coeff1) in p1.ops
        for (op2,coeff2) in p2.ops
            out.ops[op1 ⊗ op2] = coeff1 * coeff2 
        end
    end
    return out 
end
const ⊗ = otimes
