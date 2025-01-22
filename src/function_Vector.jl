

"""
    Base.Vector(k::SparseKetBasis{N,T}) where {N,T}

Create dense vector representation in standard basis 
"""
function Base.Vector(k::SparseKetBasis{N,T}) where {N,T}
    vec = zeros(T,Int128(2)^N)
    for (ket, coeff) in k
        vec[ket.v+1] = coeff
    end
    return vec 
end


"""
    Base.Vector(k::Union{KetBitString{N}, Dyad{N}, FixedPhasePauli{N}) where N

Create dense vector representation in standard basis 
"""
function Base.Vector(k::KetBitString{N}; T=Int64) where N
    vec = zeros(T,Int128(2)^N)
    vec[index(k)] = T(1) 
    return vec 
end


"""
    Base.Vector(k::Union{KetBitString{N}, Dyad{N}, FixedPhasePauli{N}) where N

Create dense vector representation in standard basis 
"""
function Base.Vector(k::Union{Dyad{N}, FixedPhasePauli{N}}; T=Int64) where N
    return vec(Matrix(k))
end