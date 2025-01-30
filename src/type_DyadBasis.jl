"""
A basis for `Dyad`'s, which is not closed under multiplication. 
Since the product of two arbitrary dyad's don't generally create another dyad (e.g., while |i><j| * |j><l| = |i><l|, most products create scalars: |i><j| * |k><l| = 0).
As such the product of two `DyadBasis` objects is not a `DyadBasis` object, but a `Dyad`, which contains a scalar factor. 
This type is primarily used to provide a basis for linear combinations of `Dyad`'s, e.g., `DyadSum`'s.
"""
struct DyadBasis{N}  
    ket::Ket{N}
    bra::Bra{N} 
end


DyadBasis(d::DyadBasis) = d

"""
    DyadBasis(ket::Vector{T}, bra::Vector{T}) where T<:Union{Bool, Integer}

TBW
"""
function DyadBasis(ket::Vector{T}, bra::Vector{T}) where T<:Union{Bool, Integer}
    @assert length(ket) == length(bra)
    N = length(ket) 
    return DyadBasis{N}(Ket(ket),Bra(bra))
end

"""
    DyadBasis(N::Integer, k::Integer, b::Integer)

TBW
"""
function DyadBasis(N::Integer, k::Integer, b::Integer)
    return DyadBasis{N}(Ket{N}(Int128(k)), Bra{N}(Int128(b)))
end

Base.rand(::Type{DyadBasis{N}}) where N = DyadBasis{N}(rand(Ket{N}), rand(Bra{N}))
          
Base.adjoint(d::DyadBasis{N}) where N = DyadBasis{N}(adjoint(d.bra), adjoint(d.ket))

@inline LinearAlgebra.ishermitian(d::DyadBasis) = d.ket.v == d.bra.v


function Base.iterate(::Type{DyadBasis{N}}, state = 1) where N
    state > 4^N && return
    next = CartesianIndices((2^N,2^N))[state]
    return DyadBasis{N}(Ket{N}(next[1]-1), Bra{N}(next[2]-1)), state+1 
end

Base.display(d::DyadBasis{N}) where N = println(string(d))
function Base.string(d::DyadBasis{N}) where N
    return string(d.ket)*string(d.bra)
end

coeff(d::DyadBasis) = true 


function otimes(p1::DyadBasis{N}, p2::DyadBasis{M}) where {N,M} 
    DyadBasis{N+M}(p1.ket ⊗ p2.ket, p1.bra ⊗ p2.bra)
end