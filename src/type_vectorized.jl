"""
    ps::PauliSum{N}

A vectorized form of `PauliSum`. Essentially a Liouville space form.
"""
struct VectorizedPauliSum{N} <: AbstractArray{ComplexF64,1} 
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


struct SuperOperator{T,N} <: AbstractMatrix{T}
    matvec
    dim::Int
    sym::Bool
end

Base.size(lop::SuperOperator{T}) where {T} = return (lop.dim,lop.dim)
Base.:(*)(lop::SuperOperator{T}, v::AbstractVector{T}) where {T} = return lop.matvec(v)
# Base.:(*)(lop::SuperOperator{T}, v::AbstractMatrix{T}) where {T} = return lop.matvec(v)
issymmetric(lop::SuperOperator{T}) where {T} = return lop.sym
    
function Base.display(L::SuperOperator{T,N}) where {T,N}
    @printf("SuperOperator: dim = %5i N = %2i\n", L.dim, N)
end
Base.show(L::SuperOperator) = display(L)

function vectorized_commutator(H::PauliSum{N}; T=ComplexF64) where N

    function mymatvec(v::VectorizedPauliSum{N})
        return VectorizedPauliSum(commutator(H, v.ps))
    end

    return SuperOperator{T, N}(mymatvec, 2^N, false)
end


function vectorized_rmul(H::PauliSum{N}; T=ComplexF64) where N

    function mymatvec(v::VectorizedPauliSum{N})
        return VectorizedPauliSum(v.ps * H)
    end

    return SuperOperator{T, N}(mymatvec, 2^N, false)
end


function vectorized_lmul(H::PauliSum{N}; T=ComplexF64) where N

    function mymatvec(v::VectorizedPauliSum{N})
        return VectorizedPauliSum(H * v.ps)
    end

    return SuperOperator{T, N}(mymatvec, 2^N, false)
end

