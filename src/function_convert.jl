
function Base.convert(::Type{Pauli{N}}, p::PauliPF{N}) where {N}
    return Pauli{N}(phase(p), p.z, p.x)
end

function Pauli(p::PauliPF{N}) where {N}
    return Pauli{N}(phase(p), p.z, p.x)
end

function Base.convert(::Type{ScaledPauli{T,N}}, p::Pauli{N}) where {T,N}
    return ScaledPauli{T,N}(get_phase(p), PauliPF{N}(p))
end

function Base.convert(::Type{ScaledPauli{T,N}}, p::PauliPF{N}) where {T,N}
    return ScaledPauli{T,N}(T(1), p)
end

phasefree(p::PauliPF) = p 
phasefree(p::Pauli{N}) where N = PauliPF{N}(p.z,p.x) 
phasefree(p::ScaledPauli) = phasefree(p.pauli) 