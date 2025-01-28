
"""
    PauliBasis(p::Pauli{N}) where N

Return the `PauliBasis` with the same operator part as `p`
"""
function PauliBasis(p::Pauli{N}) where N
    return PauliBasis{N}(p.z, p.x)
end
PauliSum(p::Pauli{N}; T=ComplexF64) where {N} = PauliSum{N, T}(PauliBasis(p)=>1im^global_phase(p))
PauliSum(p::PauliBasis{N}; T=ComplexF64) where {N} = PauliSum{N, T}(PauliBasis(p)=>1im^global_phase(p))

function DyadBasis(d::Dyad{N}) where N
    return DyadBasis{N}(d.ket, d.bra)
end

function Dyad(d::DyadBasis{N}) where N
    return Dyad{N}(true, d.ket, d.bra)
end