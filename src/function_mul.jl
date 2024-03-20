using InteractiveUtils

"""
    mult!(out::Matrix{T}, p::Pauli{N}, in::Matrix{T}) where {T,N}

Multiply dense vectors in `in` with `p`, storing the result in `out`. 

I think this should be much faster if we were to store vectors as rows, 
so that summing over states acts on contiguous data 
"""
function LinearAlgebra.mul!(out::Matrix{T}, p::AbstractPauli{N}, in::Matrix{T}) where {T,N}
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
    LinearAlgebra.mul!(out::Matrix{T}, p::PauliSum{N}, in::Matrix{T}) where {T,N}

Multiply dense vectors in `in` with `p`, storing the result in `out`. 

I think this should be much faster if we were to store vectors as rows, 
so that summing over states acts on contiguous data 
"""
function LinearAlgebra.mul!(out::Matrix{T}, p::PauliSum{N}, in::Matrix{T}) where {T,N}
    # fill!(out, T(0))
    ndim = size(in,1)
    nvec = size(in,2)
  
    # if nvec == 1
    #     mul!(@view out[:,1], p, in[:,1])
    #     return
    # end
    
    # Check dimensions 
    size(in) == size(out) || throw(DimensionMismatch)
    ndim == 2^N || throw(DimensionMismatch)

    if nvec == 1
        for (op,coeff) in p.ops
            # mul!(out, op, in, coeff, 1.0)
            @inbounds @simd for i in 0:ndim-1                           
                (phase, j) = op * KetBitString{N}(i)
                tmp = phase * coeff 
                k = 1
                out[j.v+1] += tmp * in[i+1]
            end
        end
    elseif nvec > 1
        for (op,coeff) in p.ops
            # mul!(out, op, in, coeff, 1.0)
            @inbounds @simd for i in 0:ndim-1                           
                (phase, j) = op * KetBitString{N}(i)
                tmp = phase * coeff 
                for k in 1:nvec
                    # out[j.v+1 + k*ndim] += tmp * in[i+1, k]
                    out[j.v+1 + (k-1)*ndim] += tmp * in[i+1 + (k-1)*ndim]
                end
            end
        end
    end
end


# function mul2!(C::Matrix{T}, A::PauliSum{N}, B::Matrix{T}) where {T,N}
#     ndim = size(B,1)
#     nvec = size(B,2)
  
#     _B = zeros(T, size(B'))
#     _C = zeros(T, size(C'))
#     _B .= transpose(B)
#     _C .= transpose(C)
   
#     # Check dimensions 
#     size(B) == size(C) || throw(DimensionMismatch)
#     ndim == 2^N || throw(DimensionMismatch)

#     for (op,coeff) in A.ops
#         # mul!(out, op, in, coeff, 1.0)
#         for i in 0:ndim-1                           
#             (phase, j) = op * KetBitString{N}(i)
#             tmp = phase * coeff 
#             # _C[:, j.v+1] .+= tmp .* _B[:, i+1]
#             @inbounds @simd for k in 1:nvec
#                 _C[k, j.v+1] += tmp * _B[k, i+1]
#             end
#         end
#     end
#     C .= transpose(_C)
# end

function LinearAlgebra.mul!(C::Vector{T}, A::PauliSum{N}, B::Vector{T}) where {T,N}
    ndim = size(B,1)
  
    # Check dimensions 
    size(B) == size(C) || throw(DimensionMismatch)
    ndim == 2^N || throw(DimensionMismatch)

    for (op,coeff) in A.ops
        @inbounds @simd for i in 0:ndim-1                           
            (phase, j) = op * KetBitString{N}(i)
            C[j.v+1] += phase * coeff * B[i+1]
        end
    end
end


function LinearAlgebra.mul!(C::Vector{T}, A::ScaledPauliVector{N}, B::Vector{T}) where {T,N}
    ndim = size(B,1)
  
    # Check dimensions 
    size(B) == size(C) || throw(DimensionMismatch)
    ndim == 2^N || throw(DimensionMismatch)

    for op in A
        @inbounds @simd for i in 0:ndim-1                           
            (coeff, j) = op * KetBitString{N}(i)
            C[j.v+1] += coeff * B[i+1]
        end
    end
end


"""
    LinearAlgebra.mul!(C::Matrix{T}, A::Pauli{N}, B::Matrix{T}, α, β) where {T,N}

    mul!(C, A, B, α, β) -> C

ABα+Cβ. The result is stored in C by overwriting it. Note that C must not be aliased with either A or B.
"""
function LinearAlgebra.mul!(C::Matrix{T}, A::AbstractPauli{N}, B::Matrix{T}, α, β) where {T,N}
    # scale!(C, β)
    # for i in 1:size(C,1)
    #     # C[i,:] .= C[i,:] .* β
    # end

    C .*= β
    ndim = size(B,1)
    nvec = size(B,2)
    
    # Check dimensions 
    size(B) == size(C) || throw(DimensionMismatch)
    ndim == 2^N || throw(DimensionMismatch)

    # loop over states and do multiplication
    # need to go from 0:N-1 because we convert to binary
    for i in 0:ndim-1                           
        (coeff, j) = A * KetBitString{N}(i)
        # @views C[j.v+1,:] .+= (coeff * α) .* B[i+1,:] 
        tmp = coeff * α 
        @simd for k in 1:nvec
            @inbounds C[j.v+1, k] += tmp * B[i+1, k]
        end
    end
end


"""
    Base.:*(p::PauliSum{N}, in::Array{T}) where {T,N}

TBW
"""
function Base.:*(p::PauliSum{N}, in::Array{T}) where {T,N}
    out = zeros(T, size(in))
    mul!(out, p, in)
    return out
end
function Base.:*(p::ScaledPauliVector{N}, in::Array{T}) where {T,N}
    out = zeros(T, size(in))
    mul!(out, p, in)
    return out
end


"""
    Base.:*(p::Pauli{N}, in::Array{T}) where {T,N}

TBW
"""
function Base.:*(p::AbstractPauli{N}, in::Array{T}) where {T<:Complex,N}
    out = zeros(T, size(in))
    mul!(out, p, in)
    return out
end


"""
    Base.:*(p1::Pauli{N}, p2::Pauli{N}) where {N}

Multiply two `Pauli`'s together
"""
function Base.:*(p1::Pauli{N}, p2::Pauli{N}) where {N}
    θ = (p1.θ + p2.θ + phase(p1.pauli, p2.pauli)) % 4
    return Pauli{N}(θ, p1.pauli*p2.pauli)
end


function Base.:*(p1::Pauli{N}, p2::FixedPhasePauli{N}) where {N}
    # x = p1.x ⊻ p2.x
    # z = p1.z ⊻ p2.z
    # θ = (p1.θ + phase(p2)) % 4
    # θ += (2*count_ones(p1.x & p2.z)) % 4
    # return Pauli{N}(θ, z, x)
    θ = (p1.θ + phase(p1.pauli, p2)) % 4
    return Pauli{N}(θ, p1.pauli*p2)
end

function Base.:*(p1::FixedPhasePauli{N}, p2::Pauli{N}) where {N}
    # x = p1.x ⊻ p2.x
    # z = p1.z ⊻ p2.z
    # θ = (phase(p1) + p2.θ) % 4
    # θ += (2*count_ones(p1.x & p2.z)) % 4
    # return Pauli{N}(θ, z, x)
    θ = (p2.θ + phase(p1, p2.pauli)) % 4
    return Pauli{N}(θ,p1*p2.pauli)
end

Base.:*(p1::FixedPhasePauli{N}, p2::ScaledPauli{N}) where N = ScaledPauli{N}(p1 * p2.coeff, p1 * p2.pauli)
Base.:*(p1::ScaledPauli{N}, p2::FixedPhasePauli{N}) where N = ScaledPauli{N}(p1.coeff, p1.pauli * p2)
Base.:*(p1::ScaledPauli{N}, p2::ScaledPauli{N}) where N = ScaledPauli{N}(p1.coeff*p2.coeff * get_phase(p1.pauli, p2.pauli), p1.pauli*p2.pauli)
Base.:*(p1::ScaledPauli{N}, p2::Pauli{N}) where N = ScaledPauli{N}(p1.coeff * get_phase(p1.pauli, phasefree(p2)), p1.pauli*p2)
Base.:*(p1::Pauli{N}, p2::ScaledPauli{N}) where N = ScaledPauli{N}(p2.coeff * get_phase(p2.pauli, phasefree(p1)), p1*p2.pauli)


"""
    Base.:*(p::Pauli{N}, c::Number) where {N}

Multiply a `Pauli` with a number. This returns a `PauliSum` 
"""
Base.:*(p::Pauli{N}, c::T) where {N,T<:Number}             = ScaledPauli{N}(c*get_phase(p), p.pauli)
Base.:*(p::FixedPhasePauli{N}, c::T) where {N,T<:Number}   = ScaledPauli{N}(c, p)
Base.:*(p::ScaledPauli{N}, c::T) where {N,T<:Number}       = ScaledPauli{N}(p.coeff*c, p.pauli) 
Base.:*(c::Number, p::Pauli) = p*c
Base.:*(c::Number, p::FixedPhasePauli) = p*c
Base.:*(c::Number, p::ScaledPauli) = p*c




# function Base.:*(p::ScaledPauli{N}, a::Number) where {T,N}
#     return ScaledPauli{N}(p.coeff*a, p.pauli)
# end

# Base.:*(a::Number, p::ScaledPauli{N}) where {T,N} = p*a


"""
    Base.:*(p::Pauli{N}, KetBitString{N}) where N

TBW
"""
function Base.:*(p::Pauli{N}, ψ::KetBitString{N}) where N
    tmp = p.pauli.x ⊻ ψ.v
    sign = count_ones(p.pauli.z & tmp) % 2
    return get_phase(p)*(-1)^sign, KetBitString{N}(tmp)
end
function Base.:*(p::FixedPhasePauli{N}, ψ::KetBitString{N}) where N
    tmp = p.x ⊻ ψ.v
    # sign = count_ones(p.pauli.z & tmp) % 2
    return iseven(count_ones(p.z & tmp)) ? 1 : -1, KetBitString{N}(tmp)
end
function Base.:*(p::ScaledPauli{N}, ψ::KetBitString{N}) where N
    tmp = p.pauli.x ⊻ ψ.v
    sign = count_ones(p.pauli.z & tmp) % 2
    return p.coeff*(-1)^sign, KetBitString{N}(tmp)
end

"""
    Base.:*(p::AbstractPauli{N}, ψ::SparseKetBasis{N,T}) where {N,T}

TBW
"""
function Base.:*(p::AbstractPauli{N}, ψ::SparseKetBasis{N,T}) where {N,T}

    σ = SparseKetBasis(N, T=ComplexF64)
    for (ket,coeff) in ψ.coeffs
        coeff2, ket2 = p * ket
        sum!(σ, ket2, coeff2*coeff)
    end
    return σ
end


function Base.:*(p::PauliSum{N}, ψ::SparseKetBasis{N,T}) where {N,T}

    σ = SparseKetBasis(N, T=ComplexF64)
    for (pauli,coeff0) in p.ops
        for (ket,coeff1) in ψ.coeffs
            coeff2, ket2 = pauli * ket
            sum!(σ, ket2, coeff2*coeff1*coeff0)
        end
    end
    return σ
end

function Base.:*(v1::SparseKetBasis{N,T}, a::Number) where {N,T}
    out = deepcopy(v1)
    scale!(out,a)
    return out
end
Base.:*(a::Number, v1::SparseKetBasis{N,T}) where {N,T} = v1 * a