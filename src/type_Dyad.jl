"""
An occupation number vectors, up to 128 qubits
"""
struct Dyad{N}  
    weight::Bool # 0 or 1, needed to close multiplication 
    ket::Ket{N}
    bra::Bra{N} 
end

"""
    Dyad(ket::Vector{T}, bra::Vector{T}) where T<:Union{Bool, Integer}

TBW
"""
function Dyad(ket::Vector{T}, bra::Vector{T}) where T<:Union{Bool, Integer}
    N = length(ket) 
    return Dyad{N}(true, Ket(ket),Bra(bra))
end

"""
    Dyad(N::Integer, k::Integer, b::Integer)

TBW
"""
function Dyad(N::Integer, k::Integer, b::Integer)
    return Dyad{N}(true, Ket{N}(k), Bra{N}(b))
end

@inline coeff(d::Dyad) = d.weight 
          
Base.adjoint(d::Dyad{N}) where N = Dyad{N}(coeff(d), adjoint(d.bra), adjoint(d.ket))

function Base.:*(d1::Dyad{N}, d2::Dyad{N}) where N
    return Dyad{N}(coeff(d1) & coeff(d2) & (d1.bra.v==d2.ket.v), d1.ket, d2.bra)
end

function Base.Matrix(d::Union{Dyad{N}, DyadBasis{N}}) where N
    mat = zeros(Bool, size(d))
    mat[d.ket.v+1, d.bra.v+1] = coeff(d) 
    return mat 
end

function Base.size(d::Union{Dyad{N}, DyadBasis{N}}) where N
    return (BigInt(2)^N, BigInt(2)^N)
end

Base.rand(::Type{Dyad{N}}) where N = Dyad{N}(true, rand(Ket{N}), rand(Bra{N}))

@inline is_hermitian(d::Dyad) = d.ket.v == d.bra.v

"""
    Base.:+(p::Dyad{N}, q::Dyad{N}) where N

Add two `Dyad`'s together, return a `DyadSum`
"""
function Base.:+(p::Dyad{N}, q::Dyad{N}) where N
    if DyadBasis(p) == DyadBasis(q)
        return DyadSum{N, ComplexF64}(DyadBasis(p) => coeff(p)+coeff(q))
    else 
        return DyadSum{N, ComplexF64}(DyadBasis(p) => coeff(p), DyadBasis(q) => coeff(q))
    end
end

function Base.iterate(::Type{Dyad{N}}, state = 1) where N
    state > 4^N && return
    next = CartesianIndices((2^N,2^N))[state]
    return Dyad{N}(true, Ket{N}(next[1]-1), Bra{N}(next[2]-1)), state+1 
end

Base.display(d::Dyad{N}) where N = println(string(d))
function Base.string(d::Dyad{N}) where N
    return "|"*string(d.ket.v)*"><"*string(d.bra.v)*"|"
end

LinearAlgebra.tr(d::Union{Dyad{N}, DyadBasis{N}}) where N = coeff(d)  && (d.ket.v == d.bra.v) 
