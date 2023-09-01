"""
    coeff::T
    pauli::Pauli{N}
"""
struct ScaledPauli{T,N}
    coeff::T
    pauli::Pauli{N}
end

function Base.convert(::Type{ScaledPauli{T,N}}, p::Pauli{N}) where {T,N}
    return ScaledPauli{T,N}(T(0), p)
end

function ScaledPauli(p::Pauli{N}) where N
    return ScaledPauli{Float64, N}(1.0, p) 
end

function Base.display(sp::ScaledPauli)
    @printf(" %12.8f * %1i+%1ii %s\n", sp.coeff, real(get_phase(sp.pauli)), imag(get_phase(sp.pauli)), string(sp.pauli))
end


"""
    Base.:-(p::ScaledPauli{T,N}) where {T,N}

Negate the `p`. We could either change the coeff or the pauli. Not sure which is best, negate coeff currently
"""
function Base.:-(p::ScaledPauli{T,N}) where {T,N}
    return ScaledPauli{T,N}(-p.coeff, p.pauli)
end