abstract type AbstractState{N} end

"""
An occupation number vector, up to 128 qubits
"""
struct KetBitString{N} <: AbstractState{N} 
    v::Int128
end
struct BraBitString{N} <: AbstractState{N}
    v::Int128
end
    
function Base.size(d::KetBitString{N}) where N
    return (BigInt(2)^N, 1)
end

function Base.size(d::BraBitString{N}) where N
    return (1, BigInt(2)^N)
end

Base.adjoint(d::KetBitString{N}) where N = BraBitString{N}(d.v)
Base.adjoint(d::BraBitString{N}) where N = KetBitString{N}(d.v)


Base.rand(T::Type{KetBitString{N}}) where N = T(rand(0:Int128(2)^N-1))
Base.rand(T::Type{BraBitString{N}}) where N = T(rand(0:Int128(2)^N-1))

SparseKetBasis{N, T} = Dict{KetBitString{N}, T}
function SparseKetBasis(N; T=Float64)
    return Dict{KetBitString{N}, T}()
end



"""
    KetBitString(vec::Vector{T}) where T<:Union{Bool, Integer}

TBW
"""
function KetBitString(vec::Vector{T}) where T<:Union{Bool, Integer}
    two = Int128(2)
    v = Int128(0)

    for i in 1:length(vec)
        if vec[i] == 1
            v |= two^(i-1)
        end
    end
    return KetBitString{length(vec)}(v)
end
function BraBitString(vec::Vector{T}) where T<:Union{Bool, Integer}
    two = Int128(2)
    v = Int128(0)

    for i in 1:length(vec)
        if vec[i] == 1
            v |= two^(i-1)
        end
    end
    return BraBitString{length(vec)}(v)
end

"""
    KetBitString(N::Integer, v::Integer)

TBW
"""
function KetBitString(N::Integer, v::Integer)
    for i in N+1:128
        v &= ~(Int128(2)^(i-1))    
    end
    return KetBitString{N}(Int128(v))
end
function BraBitString(N::Integer, v::Integer)
    for i in N+1:128
        v &= ~(Int128(2)^(i-1))    
    end
    return BraBitString{N}(Int128(v))
end


"""
    Base.show(io::IO, P::Pauli{N}) where N

TBW
"""
function Base.show(io::IO, P::AbstractState{N}) where N
    print(io, string(P))
end

"""
    Base.show(io::IO, P::Pauli{N}) where N

TBW
"""
function Base.show(io::IO, v::SparseKetBasis{N,T}) where {N,T}
    for (ket,coeff) in v
        print(io, string(ket), coeff)
    end
end

function Base.string(p::KetBitString{N}) where N
    out = [0 for i in 1:128]
    for i in get_on_bits(p.v)
        out[i] = 1
    end
    return "|"*join(out[1:N])*">"
end
function Base.string(p::BraBitString{N}) where N
    out = [0 for i in 1:128]
    for i in get_on_bits(p.v)
        out[i] = 1
    end
    return "<"*join(out[1:N])*"|"
end




"""
    LinearAlgebra.dot(v1::SparseKetBasis{N,T}, v2::SparseKetBasis{N,TT}) where {N,T,TT}

TBW
"""
function LinearAlgebra.dot(v1::SparseKetBasis{N,T}, v2::SparseKetBasis{N,TT}) where {N,T,TT}
    out = 0.0
    if length(v1) < length(v2)
        for (ket,coeff) in v1
            out += adjoint(coeff) * get(v2, ket, 0.0)
        end
    else
        for (ket,coeff) in v2
            out += coeff * adjoint(get(v1, ket, 0.0))
        end
    end
    return out
end

"""
    scale!(v1::SparseKetBasis{N,T}, a::Number) where {N,T}

TBW
"""
function scale!(v1::SparseKetBasis{N,T}, a::Number) where {N,T}
    map!(x->x*a, values(v1))
end

