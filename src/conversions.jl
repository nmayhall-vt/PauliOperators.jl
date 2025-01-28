
"""
    PauliBasis(p::Pauli{N}) where N

Return the `PauliBasis` with the same operator part as `p`
"""
function PauliBasis(p::Pauli{N}) where N
    return PauliBasis{N}(p.z, p.x)
end

function DyadBasis(d::Dyad{N}) where N
    return DyadBasis{N}(d.ket, d.bra)
end

function Dyad(d::DyadBasis{N}) where N
    return Dyad{N}(true, d.ket, d.bra)
end