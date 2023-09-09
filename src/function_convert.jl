
function Base.convert(::Type{Pauli{N}}, p::FixedPhasePauli{N}) where {N}
    return Pauli{N}(phase(p), p.z, p.x)
end

function Pauli(p::FixedPhasePauli{N}) where {N}
    return Pauli{N}(0, p.z, p.x)
end

function Base.convert(::Type{ScaledPauli{N}}, p::Pauli{N}) where {T,N}
    return ScaledPauli{N}(get_phase(p), FixedPhasePauli{N}(p))
end

function Base.convert(::Type{FixedPhasePauli{N}}, p::Pauli{N}) where {T,N}
    return FixedPhasePauli{N}(p.z, p.x)
end

function Base.convert(::Type{ScaledPauli{N}}, p::FixedPhasePauli{N}) where {T,N}
    return ScaledPauli{N}(T(1), p)
end

phasefree(p::FixedPhasePauli) = p 
phasefree(p::Pauli{N}) where N = FixedPhasePauli{N}(p.z,p.x) 
phasefree(p::ScaledPauli) = phasefree(p.pauli) 