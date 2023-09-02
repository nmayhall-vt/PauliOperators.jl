Base.isless(p1::Pauli{N}, p2::Pauli{N}) where N = isless((p1.z, p1.x), (p2.z, p2.x))
Base.isless(p1::ScaledPauli{T,N}, p2::Pauli{N}) where {T,N} = isless(p1.pauli, p2)
Base.isless(p1::Pauli{N}, p2::ScaledPauli{T,N}) where {T,N} = isless(p1, p2.pauli)

"""
    Base.length(spv::ScaledPauliVector)

TBW
"""
Base.length(spv::ScaledPauliVector) = length(spv.paulis)