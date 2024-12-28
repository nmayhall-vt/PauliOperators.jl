"""
An occupation number vectors, up to 128 qubits
"""
struct Dyad{N} 
    ket::KetBitString{N}
    bra::BraBitString{N} # should this be a different type?
end
          
struct ScaledDyad{N,T}
    coeff::T 
    dyad::Dyad{N}
end
          

DyadSum{N, T} = Dict{Dyad{N}, T}
function DyadSum(N; T=Float64)
    return Dict{Dyad{N}, T}()
end

Base.adjoint(d::Dyad{N}) where N = Dyad{N}(adjoint(d.bra), adjoint(d.ket))
Base.adjoint(d::ScaledDyad{N,T}) where {N,T} = ScaledDyad{N,T}(adjoint(d.coeff), adjoint(d.dyad))

"""
    Dyad(ket::Vector{T}, bra::Vector{T}) where T<:Union{Bool, Integer}

TBW
"""
function Dyad(ket::Vector{T}, bra::Vector{T}) where T<:Union{Bool, Integer}
    @assert length(ket) == length(bra)
    N = length(ket) 
    return Dyad{N}(KetBitString(ket),BraBitString(bra))
end

"""
    Dyad(N::Integer, k::Integer, b::Integer)

TBW
"""
function Dyad(N::Integer, k::Integer, b::Integer)
    return Dyad{N}(KetBitString{N}(Int128(k)), BraBitString{N}(Int128(b)))
end

"""
    Dyad(N::Integer, k::Integer, b::Integer)

TBW
"""
function ScaledDyad(N::Integer, c::T, k::Integer, b::Integer) where T
    return ScaledDyad{N,T}(c, Dyad(N, k,b))
end

function Base.rand(T::Type{Dyad{N}}) where N
    return Dyad{N}(rand(KetBitString{N}), rand(BraBitString{N}))
end

function Base.rand(T::Type{ScaledDyad{N,TT}}) where {N,TT}
    return ScaledDyad{N,TT}(rand(TT), rand(Dyad{N}))
end

function Base.rand(T::Type{ScaledDyad{N}}) where N
    return ScaledDyad{N,ComplexF64}(rand(ComplexF64), rand(Dyad{N}))
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
    return "|"*string(d.ket.v)*"><"*string(d.bra.v)*"|"
    # return "|"*string(d.ket)*"><"*string(d.bra)*"|"
    # return "|"*join(get_on_bits(d.ket.v))*"><"*join(get_on_bits(d.bra.v))*"|"
end
function Base.string(d::ScaledDyad)
    return string(d.coeff)*string(d.dyad)
end


# """
#     Base.Vector(k::KetBitString{N}) where N

# TBW
# """
# function Base.Vector(k::KetBitString{N}) where N
#     vec = zeros(Int8,Int128(2)^N)
#     vec[k.v+1] = 1
#     return vec 
# end
# """
#     Base.Vector(k::KetBitString{N}) where N

# TBW
# """
# function Base.Vector(k::SparseKetBasis{N,T}) where {N,T}
#     vec = zeros(T,Int128(2)^N)
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


function Base.Matrix(d::Dyad{N}) where N
    mat = zeros(Bool, size(d))
    mat[d.ket.v+1, d.bra.v+1] = true
    return mat 
end

function Base.size(d::Dyad{N}) where N
    return (BigInt(2)^N, BigInt(2)^N)
end
Base.size(d::ScaledDyad) = size(d.dyad) 

function Base.Matrix(d::ScaledDyad{N,T}) where {N,T}
    mat = zeros(T, size(d))
    mat[d.dyad.ket.v+1, d.dyad.bra.v+1] = d.coeff 
    return mat 
end