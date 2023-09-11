"""
    commute(p1,p2)

Check if `p1` and `p2` commute, return true or false
"""
@inline commute(p1::FixedPhasePauli, p2::FixedPhasePauli) = iseven(count_ones(p1.x & p2.z) - count_ones(p1.z & p2.x)) 
@inline commute(p1::Pauli{N}, p2::Pauli{N}) where {N} = commute(p1.pauli, p2.pauli)
@inline commute(p1::ScaledPauli{N}, p2::ScaledPauli{N}) where {N} = commute(p1.pauli, p2.pauli)
@inline commute(p1::Pauli{N}, p2::ScaledPauli{N}) where {N} = commute(p1, p2.pauli)
@inline commute(p1::ScaledPauli{N}, p2::Pauli{N}) where {N} = commute(p1.pauli, p2)