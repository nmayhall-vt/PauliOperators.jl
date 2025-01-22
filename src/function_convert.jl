
function Base.convert(::Type{Pauli{N}}, p::FixedPhasePauli{N}) where {N}
    return Pauli{N}(0, p)
end

function Base.convert(::Type{ScaledPauli{N}}, p::Pauli{N}) where {N}
    return ScaledPauli{N}(get_phase(p), p.pauli)
end

function Base.convert(::Type{ScaledPauli{N}}, p::FixedPhasePauli{N}) where {N}
    return ScaledPauli{N}(ComplexF64(1), p)
end

phasefree(p::FixedPhasePauli) = p 
phasefree(p::Pauli) = p.pauli 
phasefree(p::ScaledPauli) = p.pauli