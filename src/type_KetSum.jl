
Base.adjoint(d::KetSum{N,T}) where {N,T} = Adjoint(d)
Base.parent(d::Adjoint{<:Any, <:KetSum}) = d.parent


"""
    Base.show(io::IO, v::KetSum{N,T}) where {N,T}

TBW
"""
function Base.show(io::IO, v::KetSum{N,T}) where {N,T}
    for (ket,coeff) in v
        print(io, string(ket), coeff)
    end
end

"""
    LinearAlgebra.dot(v1::KetSum{N,T}, v2::KetSum{N,TT}) where {N,T,TT}

TBW
"""
function LinearAlgebra.dot(v1::KetSum{N,T}, v2::KetSum{N,TT}) where {N,T,TT}
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
    scale!(v1::KetSum{N,T}, a::Number) where {N,T}

TBW
"""
function scale!(v1::KetSum{N,T}, a::Number) where {N,T}
    map!(x->x*a, values(v1))
end


"""
    Base.Vector(k::KetSum{N,T}) where {N,T}

TBW
"""
function Base.Vector(k::KetSum{N,T}) where {N,T}
    vec = zeros(T,Int128(2)^N)
    for (k,coeff) in k
        vec[index(k)] = T(coeff) 
    end
    return vec 
end


function Base.Vector(k::Adjoint{<:Any, KetSum{N,T}}) where {N,T}
    vec = zeros(T,Int128(2)^N)
    for (k,coeff) in k.parent
        vec[index(k)] = T(coeff') 
    end
    return vec 
end