"""
An occupation number vectors, up to 128 qubits
"""
struct Dyad{N}  
    s::ComplexF64  
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

@inline coeff(d::Dyad) = d.s 
          
Base.adjoint(d::Dyad{N}) where N = Dyad{N}(coeff(d), adjoint(d.bra), adjoint(d.ket))

Base.:*(d1::Dyad{N}, a::Number) where N = Dyad{N}(coeff(d1)*a, d1.ket, d1.bra)
Base.:*(a::Number, d1::Dyad{N}) where N = d1 * a 

function Base.Matrix(d::Union{Dyad{N}, DyadBasis{N}}) where N
    mat = zeros(ComplexF64, size(d))
    mat[d.ket.v+1, d.bra.v+1] = coeff(d) 
    return mat 
end

function Base.size(d::Union{Dyad{N}, DyadBasis{N}}) where N
    return (BigInt(2)^N, BigInt(2)^N)
end

Base.rand(::Type{Dyad{N}}) where N = Dyad{N}(true, rand(Ket{N}), rand(Bra{N}))

@inline LinearAlgebra.ishermitian(d::Dyad) = d.ket.v == d.bra.v

function Base.:-(ps1::Dyad{N}) where {N}
    return Dyad{N}(-ps1.s, ps1.ket, ps1.bra) 
end

function Base.iterate(::Type{Dyad{N}}, state = 1) where N
    state > 4^N && return
    next = CartesianIndices((2^N,2^N))[state]
    return Dyad{N}(true, Ket{N}(next[1]-1), Bra{N}(next[2]-1)), state+1 
end

Base.display(d::Dyad{N}) where N = println(string(d))
function Base.string(d::Dyad{N}) where N
    return @sprintf "% .4f % .4fim | %s\n" real(d.s) imag(d.s) string(DyadBasis(d)) 
end

LinearAlgebra.tr(d::Union{Dyad{N}, DyadBasis{N}}) where N = coeff(d) * (d.ket.v == d.bra.v) 

function otimes(p1::Dyad{N}, p2::Dyad{M}) where {N,M} 
    Dyad{N+M}(coeff(p1) * coeff(p2), p1.ket ⊗ p2.ket, p1.bra ⊗ p2.bra)
end