struct StateVectorDense  
    coeffs::Vector{ComplexF64}
end

function StateVectorDense(N; nstates=1)
    return StateVectorDense(zeros(ComplexF64), 2^N, nstates)
end