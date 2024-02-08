struct StateVectorDense <: AbstractArray{ComplexF64,1}
    coeffs::Vector{ComplexF64}
end

Base.size(k::KetBitString{N}) where N = size(k.coeffs) 
Base.getindex(k::KetBitString{N}, idx) where N = k.coeffs[idx] 

function StateVectorDense(N; nstates=1)
    return StateVectorDense(zeros(ComplexF64), 2^N, nstates)
end