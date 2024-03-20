"""
An occupation number vector, up to 128 qubits
"""
struct KetBitString{N} <: AbstractArray{ComplexF64,1}
    v::Int128
end
         
Base.size(k::KetBitString{N}) where N = (2^N,)
Base.getindex(k::KetBitString{N}, idx) where N = idx == k.v


struct SparseKetBasis{N,T} <: AbstractArray{T,1}
    coeffs::OrderedDict{KetBitString{N}, T}
end

function SparseKetBasis(N::Integer; T=Float64)
    return SparseKetBasis{N,T}(OrderedDict{KetBitString{N}, T}())
end

Base.get(skb::SparseKetBasis{N,T}, ket::KetBitString{N}, i) where {N,T} = get(skb.coeffs, ket, i)
Base.size(skb::SparseKetBasis{N,T}) where {N,T} = (2^N,)
Base.getindex(k::SparseKetBasis, i) = getindex(k.coeffs, i) 
Base.setindex!(k::SparseKetBasis, i, j) = setindex!(k.coeffs, i, j) 

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

"""
    random_KetBitString(N::Integer)

TBW
"""
function random_KetBitString(N::Integer)
    return KetBitString(N, rand(1:2^N-1))
end

"""
    Base.show(io::IO, P::Pauli{N}) where N

TBW
"""
function Base.show(io::IO, P::KetBitString{N}) where N
    print(io, string(P))
end

"""
    Base.show(io::IO, P::Pauli{N}) where N

TBW
"""
function Base.show(io::IO, v::SparseKetBasis{N,T}) where {N,T}
    for (ket,coeff) in v.coeffs
        print(io, string(ket), coeff)
    end
end

"""
    Base.string(p::KetBitString{N}) where N

Display, y = iY
"""
function Base.string(p::KetBitString{N}) where N
    out = [0 for i in 1:128]
    for i in get_on_bits(p.v)
        out[i] = 1
    end
    return join(out[1:N])
end


"""
    Base.Vector(k::KetBitString{N}) where N

TBW
"""
function Base.Vector(k::KetBitString{N}) where N
    vec = zeros(Int8,2^N)
    vec[k.v+1] = 1
    return vec 
end
"""
    Base.Vector(k::KetBitString{N}) where N

TBW
"""
function Base.Vector(k::SparseKetBasis{N,T}) where {N,T}
    vec = zeros(T,2^N)
    for (ket, coeff) in k.coeffs
        vec[ket.v+1] = coeff
    end
    return vec 
end


"""
    LinearAlgebra.dot(v1::SparseKetBasis{N,T}, v2::SparseKetBasis{N,TT}) where {N,T,TT}

TBW
"""
function LinearAlgebra.dot(v1::SparseKetBasis{N,T}, v2::SparseKetBasis{N,TT}) where {N,T,TT}
    out = 0.0
    if length(v1.coeffs) < length(v2.coeffs)
        for (ket,coeff) in v1.coeffs
            out += adjoint(coeff) * get(v2, ket, 0.0)
        end
    else
        for (ket,coeff) in v2.coeffs
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
    map!(x->x*a, values(v1.coeffs))
end

