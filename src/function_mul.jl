"""
    mult!(out::Matrix{T}, p::Pauli{N}, in::Matrix{T}) where {T,N}

Multiply dense vectors in `in` with `p`, storing the result in `out`. 

I think this should be much faster if we were to store vectors as rows, 
so that summing over states acts on contiguous data 
"""
function LinearAlgebra.mul!(out::Matrix{T}, p::Pauli{N}, in::Matrix{T}) where {T,N}
    # fill!(out, T(0))
    ndim = size(in,1)
    nvec = size(in,2)
    
    # Check dimensions 
    size(in) == size(out) || throw(DimensionMismatch)
    ndim == 2^N || throw(DimensionMismatch)

    # loop over states and do multiplication

    # need to go from 0:N-1 because we convert to binary
    for i in 0:ndim-1                           
        (phase, j) = p * KetBitString{N}(i)
        # @inbounds out[:, j.v+1] .+= phase .* in[:, i+1]
        @simd for k in 1:nvec
            @inbounds out[j.v+1, k] += phase * in[i+1, k]
        end
    end
end

"""
    LinearAlgebra.mul!(C::Matrix{T}, A::Pauli{N}, B::Matrix{T}, α, β) where {T,N}

    mul!(C, A, B, α, β) -> C

ABα+Cβ. The result is stored in C by overwriting it. Note that C must not be aliased with either A or B.
"""
function LinearAlgebra.mul!(out::Matrix{T}, p::Pauli{N}, in::Matrix{T}, α, β) where {T,N}
    out .= out .* β
    ndim = size(in,1)
    
    # Check dimensions 
    size(in) == size(out) || throw(DimensionMismatch)
    ndim == 2^N || throw(DimensionMismatch)

    # loop over states and do multiplication
    # need to go from 0:N-1 because we convert to binary
    for i in 0:ndim-1                           
        (phase, j) = p * KetBitString{N}(i)
        out[j.v+1,:] .+= phase*in[i+1,:] .* α
    end
end


"""
    Base.:*(p::Pauli{N}, in::Array{T}) where {T,N}

TBW
"""
function Base.:*(p::Pauli{N}, in::Array{T}) where {T<:Complex,N}
    out = zeros(T, size(in))
    mul!(out, p, in)
    return out
end


"""
    Base.:*(p1::Pauli{N}, p2::Pauli{N}) where {N}

Multiply two `Pauli`'s together
"""
function Base.:*(p1::Pauli{N}, p2::Pauli{N}) where {N}
    x = p1.x ⊻ p2.x
    z = p1.z ⊻ p2.z
    θ = (p1.θ + p2.θ ) % 4
    θ = (θ + 2*count_ones(p1.x & p2.z)) % 4
    return Pauli{N}(θ,z,x)
end

"""
    Base.:*(p::Pauli{N}, c::Number) where {N}

Multiply a `Pauli` with a number. This returns a `PauliSum` 
"""
function Base.:*(p::Pauli{N}, c::Number) where {N}
    ps = PauliSum(N)
    sum!(ps,p)
    mul!(ps, c)
    return ps 
end
Base.:*(c::Number, p::Pauli) = p*c


"""
    Base.:*(p::Pauli{N}, KetBitString{N}) where N

TBW
"""
function Base.:*(p::Pauli{N}, ψ::KetBitString{N}) where N
    tmp = p.x ⊻ ψ.v
    sign = count_ones(p.z & tmp) % 2
    return get_phase(p)*(-1)^sign, KetBitString{N}(tmp)
end

function Base.:*(p1::ScaledPauli{T,N}, p2::ScaledPauli{T,N}) where {T,N}
    return ScaledPauli{T,N}(p1.coeff*p2.coeff, p1.pauli*p2.pauli)
end

function Base.:*(p1::ScaledPauli{T,N}, p2::Pauli{N}) where {T,N}
    return ScaledPauli{T,N}(p1.coeff, p1.pauli*p2)
end

function Base.:*(p1::Pauli{N}, p2::ScaledPauli{T,N}) where {T,N}
    return ScaledPauli{T,N}(p2.coeff, p1*p2.pauli)
end

function Base.:*(p::ScaledPauli{T,N}, a::Number) where {T,N}
    return ScaledPauli{T,N}(p.coeff*a, p.pauli)
end

Base.:*(a::Number, p::ScaledPauli{T,N}) where {T,N} = p*a

"""
    Base.:*(ps::PauliSum{N}, v::Array{T}) where {T,N}

TBW
"""
function Base.:*(ps::PauliSum{N}, v::Array{T}) where {T,N}
    σ = zeros(T,size(v))
    ndim = size(v,1)
    ndim == 2^N || throw(DimensionMismatch)
    for (op,coeff) in ps
        mul!(σ, op, v, coeff, 1.0)
    end
    return σ 
end

"""
    Base.:*(spv::Vector{ScaledPauli{T,N}}, v::Array{T}) where {T,N}

TBW
"""
function Base.:*(spv::Vector{ScaledPauli{T,N}}, v::Array{T}) where {T,N}
    σ = zeros(T,size(v))
    ndim = size(v,1)
    ndim == 2^N || throw(DimensionMismatch)
    for sp in spv 
        mul!(σ, sp.pauli, v, sp.coeff, 1.0)
    end
    return σ 
end
