"""
An occupation number vectors, up to 128 qubits
"""
struct Dyad{N} 
    ket::KetBitString{N}
    bra::KetBitString{N} # should this be a different type?
end
          
struct ScaledDyad{N,T}
    coeff::T 
    dyad::Dyad{N}
end
          

DyadSum{N, T} = Dict{Dyad{N}, T}
function DyadSum(N; T=Float64)
    return Dict{Dyad{N}, T}()
end

Base.adjoint(d::Dyad{N}) where N = Dyad{N}(d.bra,d.ket)
Base.adjoint(d::ScaledDyad{N,T}) where {N,T} = ScaledDyad{N,T}(adjoint(d.coeff)adjoint(d.dyad))

"""
    Dyad(ket::Vector{T}, bra::Vector{T}) where T<:Union{Bool, Integer}

TBW
"""
function Dyad(ket::Vector{T}, bra::Vector{T}) where T<:Union{Bool, Integer}
    @assert length(ket) == length(bra)
    N = length(ket) 
    return Dyad{N}(KetBitString(ket),KetBitString(bra))
end

"""
    Dyad(N::Integer, k::Integer, b::Integer)

TBW
"""
function Dyad(N::Integer, k::Integer, b::Integer)
    return Dyad{N}(KetBitString{N}(Int128(k)), KetBitString{N}(Int128(b)))
end

"""
    random_Dyad(N::Integer)

TBW
"""
function random_Dyad(N::Integer)
    return Dyad{N}(KetBitString(N, rand(1:2^N-1)), KetBitString(N, rand(1:2^N-1)))
end

"""
    Base.show(io::IO, P::Dyad{N}) where N

TBW
"""
function Base.show(io::IO, P::Dyad{N}) where N
    print(io, string(P))
end
Base.show(io::IO, d::ScaledDyad) = print(io,string(d))

"""
    Base.show(io::IO, v::DyadSum{N,T}) where {N,T}

TBW
"""
function Base.show(io::IO, v::DyadSum{N,T}) where {N,T}
    for (ket,coeff) in v
        print(io, string(ket), coeff)
    end
end

function Base.string(d::Dyad{N}) where N
    return "|"*join(get_on_bits(d.ket.v))*"><"*join(get_on_bits(d.bra.v))*"|"
end
function Base.string(d::ScaledDyad)
    return string(d.coeff)*string(d.dyad)
end


# """
#     Base.Vector(k::KetBitString{N}) where N

# TBW
# """
# function Base.Vector(k::KetBitString{N}) where N
#     vec = zeros(Int8,2^N)
#     vec[k.v+1] = 1
#     return vec 
# end
# """
#     Base.Vector(k::KetBitString{N}) where N

# TBW
# """
# function Base.Vector(k::SparseKetBasis{N,T}) where {N,T}
#     vec = zeros(T,2^N)
#     for (ket, coeff) in k
#         vec[ket.v+1] = coeff
#     end
#     return vec 
# end


# """
#     LinearAlgebra.dot(v1::SparseKetBasis{N,T}, v2::SparseKetBasis{N,TT}) where {N,T,TT}

# TBW
# """
# function LinearAlgebra.dot(v1::SparseKetBasis{N,T}, v2::SparseKetBasis{N,TT}) where {N,T,TT}
#     out = 0.0
#     if length(v1) < length(v2)
#         for (ket,coeff) in v1
#             out += adjoint(coeff) * get(v2, ket, 0.0)
#         end
#     else
#         for (ket,coeff) in v2
#             out += coeff * adjoint(get(v1, ket, 0.0))
#         end
#     end
#     return out
# end

# """
#     scale!(v1::SparseKetBasis{N,T}, a::Number) where {N,T}

# TBW
# """
# function scale!(v1::SparseKetBasis{N,T}, a::Number) where {N,T}
#     map!(x->x*a, values(v1))
# end

