
function Base.convert(::Type{Pauli{N}}, p::FixedPhasePauli{N}) where {N}
    return Pauli{N}(phase(p), p.pauli)
end

function Pauli(p::FixedPhasePauli{N}) where {N}
    return Pauli{N}(0, p)
end

function Base.convert(::Type{ScaledPauli{N}}, p::Pauli{N}) where {T,N}
    return ScaledPauli{N}(get_phase(p), FixedPhasePauli{N}(p))
end

function Base.convert(::Type{FixedPhasePauli{N}}, p::Pauli{N}) where {T,N}
    return p.pauli
end

function Base.convert(::Type{ScaledPauli{N}}, p::FixedPhasePauli{N}) where {T,N}
    return ScaledPauli{N}(T(1), p)
end

phasefree(p::FixedPhasePauli) = p 
phasefree(p::Pauli) = p.pauli 
phasefree(p::ScaledPauli) = p.pauli