
function Pauli(p::FixedPhasePauli{N}) where {N}
    return Pauli{N}(0, p)
end
function ScaledPauli(p::FixedPhasePauli{N}) where {N}
    return ScaledPauli{N}(1, p)
end
function ScaledPauli(p::Pauli{N}) where N
    return ScaledPauli{N}(get_phase(p), p.pauli) 
end
function PauliSum(o::FixedPhasePauli{N}) where N
    return PauliSum{N}(OrderedDict(o => 1.0))
end
function PauliSum(o::ScaledPauli{N}) where N
    return PauliSum{N}(OrderedDict(o.pauli => o.coeff))
end
function PauliSum(o::Pauli{N}) where N
    return PauliSum{N}(OrderedDict(o.pauli => get_phase(o)))
end

function Base.convert(::Type{Pauli{N}}, p::FixedPhasePauli{N}) where {N}
    return Pauli{N}(0, p.pauli)
end

function Base.convert(::Type{ScaledPauli{N}}, p::FixedPhasePauli{N}) where {N}
    return ScaledPauli{N}(1.0, p.pauli)
end

function Base.convert(::Type{ScaledPauli{N}}, p::Pauli{N}) where {T,N}
    return ScaledPauli{N}(get_phase(p), FixedPhasePauli{N}(p))
end


phasefree(p::FixedPhasePauli) = p 
phasefree(p::Pauli) = p.pauli 
phasefree(p::ScaledPauli) = p.pauli