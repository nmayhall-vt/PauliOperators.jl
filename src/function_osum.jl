

"""
    osum(p1::Pauli{N}, p2::Pauli{M}) where {N,M}

Returns the direct sum of two Paulis
"""
function osum(p1::Pauli{N}, p2::Pauli{M}) where {N,M}
    return p1 ⊗ Pauli(M) + Pauli(N) ⊗ p2 
end

"""
    osum(p1::PauliSum{N}, p2::Pauli{M}) where {N,M}

Returns the direct sum of a PauliSum and a Pauli
"""
function osum(p1::PauliSum{N}, p2::Pauli{M}) where {N,M}
    return p1 ⊗ Pauli(M) + Pauli(N) ⊗ p2 
end

"""
    osum(p1::Pauli{N}, p2::PauliSum{M}) where {N,M}

Returns the direct sum of a PauliSum and a Pauli
"""
function osum(p1::Pauli{N}, p2::PauliSum{M}) where {N,M}
    return p1 ⊗ Pauli(M) + Pauli(N) ⊗ p2 
end

"""
    osum(p1::Pauli{N}, p2::Pauli{M}) where {N,M}

Returns the direct sum of two PauliSums
"""
function osum(p1::PauliSum{N}, p2::PauliSum{M}) where {N,M}
    return p1 ⊗ Pauli(M) + Pauli(N) ⊗ p2 
end

const ⊕ = osum