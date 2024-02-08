"""
    ps::PauliSum{N}

A vectorized form of `PauliSum`. Essentially a Liouville space form.
"""
struct VectorizedPauliSum{N} <: AbstractArray{ComplexF64,1} 
    ps::PauliSum{N}
end

"""
    ps::PauliSum{N}

A vectorized form of right multiply superoperator that takes `A`` and returns `A*ps` 
"""
struct VectorizedRMult{N} <: AbstractArray{ComplexF64,2} 
    ps::PauliSum{N}
end

"""
    ps::PauliSum{N}

A vectorized form of left multiply superoperator that takes `A`` and returns `ps*A` 
"""
struct VectorizedLMult{N} <: AbstractArray{ComplexF64,2} 
    ps::PauliSum{N}
end

struct VectorizedConjugate{N} <: AbstractArray{ComplexF64,2} 
    ps::PauliSum{N}
end

struct VectorizedCommutator{N} <: AbstractArray{ComplexF64,2} 
    ps::PauliSum{N}
end

function commutator(a::PauliSum{N}, b::PauliSum{N}) where N
    out = PauliSum(N)
    for (op_a, coeff_a) in a.ops
        for (op_b, coeff_b) in b.ops
            commute(op_a, op_b) == false || continue
            out += 2*(coeff_a * op_a) * (coeff_b * op_b) 
        end
    end
    return out
end



function Base.:*(vc::VectorizedCommutator{N}, vps::VectorizedPauliSum{N}) where N
    return commutator(vc.ps, vps.ps)
end