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

# struct Adjoint{DyadSum{N,T}}
# end

DyadSum{N,T} = Dict{Dyad{N},T}
DyadSum(N::Integer; T=Float64) = return Dict{Dyad{N}, T}()
DyadSum(d::Dyad{N}; T=Float64) where N = Dict{Dyad{N}, T}(d=>T(1))
DyadSum(d::ScaledDyad{N,T}) where {N,T} = Dict{Dyad{N}, T}(d.dyad=>d.coeff)


Base.adjoint(d::Dyad{N}) where N = Dyad{N}(adjoint(d.bra), adjoint(d.ket))
Base.adjoint(d::ScaledDyad{N,T}) where {N,T} = ScaledDyad{N,T}(adjoint(d.coeff), adjoint(d.dyad))
Base.adjoint(d::DyadSum{N,T}) where {N,T} = Adjoint(d)
Base.parent(d::Adjoint{<:Any, <:DyadSum}) = d.parent

function Base.getindex(ds::Adjoint{<:Any,DyadSum{N,T}}, d::Dyad{N}) where {N,T}
    return parent(ds)[d']
end



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

function ScaledDyad(N::Integer, c::T, k::Integer, b::Integer) where T
    return ScaledDyad{N,T}(c, Dyad(N, k,b))
end
ScaledDyad(d::Dyad{N}; T=Bool) where N = ScaledDyad{N,T}(T(1), d)
ScaledDyad(c::T, d::Dyad{N}) where {T,N} = return ScaledDyad{N,T}(c,d) 

function Base.rand(T::Type{Dyad{N}}) where N
    return Dyad{N}(rand(KetBitString{N}), rand(BraBitString{N}))
end

function Base.rand(T::Type{ScaledDyad{N,TT}}) where {N,TT}
    return ScaledDyad{N,TT}(rand(TT), rand(Dyad{N}))
end

function Base.rand(T::Type{ScaledDyad{N}}) where N
    return ScaledDyad{N,ComplexF64}(rand(ComplexF64), rand(Dyad{N}))
end

function Base.:+(d1::ScaledDyad{N}, d2::ScaledDyad{N}) where {N}
    if isequal(d1.dyad, d2.dyad)
        return Dict{Dyad{N},ComplexF64}(d1.dyad=>d1.coeff+d2.coeff)
    else
        return Dict{Dyad{N},ComplexF64}(d1.dyad=>d1.coeff, d2.dyad=>d2.coeff)
    end
end
function Base.:+(a::Dyad{N}, b::ScaledDyad{N,T}) where {N,T}
    return  ScaledDyad(a, T=T) + b 
end
function Base.:+(a::ScaledDyad{N,T}, b::Dyad{N}) where {N,T}
    return  a + ScaledDyad(b, T=T)
end
function Base.:+(a::Dyad{N}, b::Dyad{N}; T=Bool) where {N}
    return  ScaledDyad(a, T=T) + ScaledDyad(b, T=T)
end
function Base.:+(a::DyadSum{N,T}, b::DyadSum{N,T}) where {N,T}
    return mergewith(+,a,b)
end
function Base.:+(a::DyadSum{N,T}, b::Adjoint{<:Any, DyadSum{N,T}}) where {N,T}
    out = deepcopy(a)
    sum!(a,b)
    return out
end
function Base.sum!(a::DyadSum{N,T}, b::DyadSum{N,T}) where {N,T}
    mergewith!(+,a,b)
end

function Base.sum!(a::DyadSum{N,T}, b::Adjoint{<:Any, DyadSum{N,T}}) where {N,T}
    for (key,val) in parent(b)
        if haskey(a, key')
            a[key'] += val'
        else
            a[key'] = val'
        end
    end
end

function Base.:+(a::Dyad{N}, b::DyadSum{N,T}) where {N,T}
    return mergewith(+,DyadSum(a,T=T),b)
end
function Base.:+(a::ScaledDyad{N}, b::DyadSum{N}) where {N}
    return mergewith(+,DyadSum(a),b)
end
function Base.:+(b::DyadSum{N},a::ScaledDyad{N}) where {N}
    return mergewith(+, b, DyadSum(a))
end

function Base.:-(d1::ScaledDyad{N}, d2::ScaledDyad{N}) where {N}
    if isequal(d1.dyad, d2.dyad)
        return Dict{Dyad{N},ComplexF64}(d1.dyad=>d1.coeff-d2.coeff)
    else
        return Dict{Dyad{N},ComplexF64}(d1.dyad=>d1.coeff, d2.dyad=>-d2.coeff)
    end
end
function Base.:-(a::Dyad{N}, b::ScaledDyad{N,T}) where {N,T}
    return  ScaledDyad(a, T=T) - b 
end
function Base.:-(a::ScaledDyad{N,T}, b::Dyad{N}) where {N,T}
    return  a - ScaledDyad(b, T=T)
end
function Base.:-(a::Dyad{N}, b::Dyad{N}; T=Bool) where {N}
    return  ScaledDyad(a, T=T) - ScaledDyad(b, T=T)
end
function Base.:-(a::DyadSum{N,T}, b::DyadSum{N,T}) where {N,T}
    return mergewith(+,a,-b)
end
Base.:-(a::Dyad{N}) where {N} = ScaledDyad{N,Bool}(false, a) 
Base.:-(a::ScaledDyad{N,T}) where {N,T} = ScaledDyad{N,T}(-a.coeff, a.dyad) 
function Base.:-(a::DyadSum{N,T}) where {N,T} 
    b = deepcopy(a)
    map!(x->-x, values(b))
    return b
end



### Multiplication 
function Base.:*(d1::Dyad{N}, d2::Dyad{N}) where {N}
    return ScaledDyad{N,Bool}(d1.bra.v==d2.ket.v, Dyad{N}(d1.ket, d2.bra))
end
function Base.:*(d1::Dyad{N}, d2::ScaledDyad{N,T}) where {N,T}
    return ScaledDyad{N,T}((d1.bra.v==d2.dyad.ket.v)*d2.coeff, Dyad{N}(d1.ket, d2.dyad.bra))
end
function Base.:*(d1::ScaledDyad{N,T}, d2::ScaledDyad{N,T}) where {N,T}
    return ScaledDyad{N,T}((d1.dyad.bra.v==d2.dyad.ket.v)*d1.coeff*d2.coeff, Dyad{N}(d1.dyad.ket, d2.dyad.bra))
end
Base.:*(a::T, d::Dyad{N}) where {N,T<:Number} = ScaledDyad{N,T}(a, d)
Base.:*(d::Dyad{N}, a::T) where {N,T<:Number} = ScaledDyad{N,T}(a, d)
Base.:*(a::T, d::ScaledDyad{N}) where {N,T<:Number} = ScaledDyad{N,T}(a*d.coeff, d.dyad)
Base.:*(d::ScaledDyad{N}, a::T) where {N,T<:Number} = ScaledDyad{N,T}(a*d.coeff, d.dyad)
Base.:*(a::T, d::DyadSum{N}) where {N,T<:Number} = replace(kv -> kv[1] => kv[2]*a, d)
Base.:*(d::DyadSum{N}, a::T) where {N,T<:Number} = replace(kv -> kv[1] => kv[2]*a, d)

"""
    Base.show(io::IO, P::Dyad{N}) where N

TBW
"""
function Base.show(io::IO, P::Dyad{N}) where N
    print(io, string(P))
end
Base.show(io::IO, d::ScaledDyad) = print(io,string(d))

# function Base.show(io::IO, v::DyadSum{N,T}) where {N,T}
#     for (key,val) in v
#         @printf(" %12.8f +%12.8fi %s\n", real(val), imag(val), key)
#     end
# end
function Base.show(io::IO, ::MIME"text/plain", a::DyadSum)
    for (key,val) in a
        show(io, "text/plain", (key,val))
        println()
    end
end


function Base.string(d::Dyad{N}) where N
    return "|"*string(d.ket.v)*"><"*string(d.bra.v)*"|"
end
function Base.string(d::ScaledDyad)
    return string(d.coeff)*string(d.dyad)
end
function Base.display(ds::DyadSum)
    for (key,val) in ds
        @printf(" %12.8f +%12.8fi %s\n", real(val), imag(val), key)
    end
end
function Base.display(ds::Adjoint{DyadSum})
    for (key,val) in ds
        @printf(" %12.8f +%12.8fi %s\n", real(val), imag(val), key)
    end
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

function Base.size(d::Dyad{N}) where N
    return (BigInt(2)^N, BigInt(2)^N)
end
Base.size(d::ScaledDyad) = size(d.dyad) 
function Base.size(d::DyadSum{N}) where N
    return (BigInt(2)^N, BigInt(2)^N)
end


function Base.Matrix(d::Dyad{N}) where N
    mat = zeros(Bool, size(d))
    mat[d.ket.v+1, d.bra.v+1] = true
    return mat 
end

function Base.Matrix(d::ScaledDyad{N,T}) where {N,T}
    mat = zeros(T, size(d))
    mat[d.dyad.ket.v+1, d.dyad.bra.v+1] = d.coeff 
    return mat 
end

function Base.Matrix(d::DyadSum{N,T}) where {N,T}
    mat = zeros(T, size(d))
    for (dyad, coeff) in d
        mat[dyad.ket.v+1, dyad.bra.v+1] = coeff 
    end
    return mat 
end