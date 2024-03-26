using LinearAlgebra
using OrderedCollections

"""
    ops::OrderedDict{Pauli{N},ComplexF64}

A collection of `Pauli`s, joined by addition.
This uses an `OrderedDict` to store them, however, the specific use cases should probably dictate the container type,
so this will probably be removed.
"""
struct PauliSum{N} <: AbstractArray{ComplexF64,2} 
    ops::OrderedDict{FixedPhasePauli{N},ComplexF64}
end

Base.size(ps::PauliSum{N}) where N = (2^N, 2^N)

"""
    PauliSum(N)

TBW
"""
function PauliSum(N::Integer)
    return PauliSum{N}(OrderedDict{FixedPhasePauli{N},ComplexF64}())
end


"""
    Base.display(ps::PauliSum)

TBW
"""
function Base.display(ps::PauliSum)
    for (key,val) in ps.ops
        @printf(" %12.8f +%12.8fi %s\n", real(val), imag(val), key)
    end
end

Base.get(ps::PauliSum{N}, p::FixedPhasePauli{N}) where N = get(ps.ops, p, zero(ComplexF64))
Base.get(ps::PauliSum{N}, p::Pauli{N}) where N = get(ps.ops, p.pauli, zero(ComplexF64))
Base.get(ps::PauliSum{N}, p::ScaledPauli{N}) where N = get(ps.ops, p.pauli, zero(ComplexF64))
Base.keys(ps::PauliSum) = keys(ps.ops)
Base.getindex(ps::PauliSum{N}, p::Pauli{N}) where N = ps.ops[p.pauli]
Base.getindex(ps::PauliSum{N}, p::FixedPhasePauli{N}) where N = ps.ops[p]
Base.setindex!(ps::PauliSum{N}, v, p::Pauli{N}) where N = ps.ops[p.pauli] = v*get_phase(p)
Base.setindex!(ps::PauliSum{N}, v, p::FixedPhasePauli{N}) where N = ps.ops[p] = v
Base.haskey(ps::PauliSum, v) = haskey(ps.ops, v)



Base.length(ps::PauliSum) = length(ps.ops)
"""
    LinearAlgebra.adjoint!(ps::PauliSum{N}) where N

TBW
"""
function LinearAlgebra.adjoint!(ps::PauliSum{N}) where N 
    for (key, val) in ps.ops
        # ps[key] = adjoint(val)  
        ps[key] = adjoint(val) * (-1)^count_ones(key.z & key.x) 
    end
end
"""
    LinearAlgebra.adjoint(ps::PauliSum{N}) where N

TBW
"""
function LinearAlgebra.adjoint(ps::PauliSum{N}) where N 
    out = deepcopy(ps)
    adjoint!(out)    
    return out
end

#
# lazy adjoint views aren't yet working. 
#
# Base.adjoint(ps::PauliSum{N}) where N = Adjoint{ComplexF64,PauliSum{N}}(ps)
# Base.size(aps::Adjoint{ComplexF64, PauliSum{N}}) where N = (length(aps.ops),)
# function Base.adjoint(ps::PauliSum)
#     return LinearAlgebra.Adjoint
# end

"""
    Base.-(p1::PauliSum{N}, p2::PauliSum{N}) where {N}

Multiply two `PauliSum`s. 
"""
function Base.:*(ps1::PauliSum{N}, ps2::PauliSum{N}) where {N}
    out = PauliSum(N)
    for (op1, coeff1) in ps1.ops 
        for (op2, coeff2) in ps2.ops
            prod = op1 * op2
            if haskey(out, prod)
                out[prod] += get_phase(op1, op2)*coeff1*coeff2
            else
                out[prod] = get_phase(op1, op2)*coeff1*coeff2
            end
            # out.ops[prod] = get(out.ops, prod) + get_phase(prod)*coeff1*coeff2 
        end
    end 
    return out
end

"""
    Base.:*(ps::PauliSum{N}, a::Number) where {N}

TBW
"""
function Base.:*(ps::PauliSum{N}, a::Number) where {N}
    out = deepcopy(ps) 
    mul!(out,a)
    return out
end
Base.:*(a::Number, ps::PauliSum{N}) where {N} = ps*a

"""
    LinearAlgebra.mul!(ps::PauliSum, a::Number)

TBW
"""
function LinearAlgebra.mul!(ps::PauliSum, a::Number)
    for (op, coeff) in ps.ops 
        ps[op] = coeff * a
    end
    return ps
end

"""
    Base.Matrix(ps::PauliSum{N}; T=ComplexF64) where N

Create a dense Matrix of type `T`
"""
function Base.Matrix(ps::PauliSum{N}; T=ComplexF64) where N
    out = zeros(T, 2^N, 2^N)
    for (op, coeff) in ps.ops
        out .+= Matrix(op) .* coeff 
    end
    return out
end

function is_hermitian(ps::PauliSum) 
    for (p,coeff) in ps.ops
        if is_hermitian(p) ⊻ (abs(imag(coeff)) < 1e-12 )
            return false
        end
    end
    return true
end

"""
    clip!(ps::PauliSum; thresh=1e-16)

Delete Pauli's with coeffs smaller than thresh
"""
function clip!(ps::PauliSum{N}; thresh=1e-16) where {N}
    filter!(p->abs(p.second) > thresh, ps.ops)
end

"""
    Base.:≈(p1::PauliSum{N}, p2::PauliSum{N}) where {N}

TBW
"""
function Base.:≈(p1::PauliSum{N}, p2::PauliSum{N}) where {N}
    for (op, coeff) in p1.ops
        if haskey(p2, op)
            coeff ≈ p2[op] || return false
        else
            return false
        end
    end 
    for (op, coeff) in p2.ops
        if haskey(p1, op)
            coeff ≈ p1[op] || return false
        else
            return false
        end
    end 
    #     get(p2, op) .≈ coeff || return false
    # end
    # for (op, coeff) in p2.ops
    #     get(p1, op) .≈ coeff || return false
    # end
    return true
end


"""
    matvec(ps::PauliSum{N}, v::Matrix) where N

TBW
"""
function matvec(ps::PauliSum{N}, v::Matrix) where N

    σ = zeros(T,size(v))

end

function LinearAlgebra.diag(ps::PauliSum{N}) where N
    out = PauliSum(N)
    for (op,coeff) in ps.ops
        if is_diagonal(op)
            out[op] = coeff
        end
    end
    return out
end

function LinearAlgebra.tr(p::PauliSum{N}) where N 
    if haskey(p, FixedPhasePauli{N}(0,0))
        return p[FixedPhasePauli{N}(0,0)] * 2^N
    else
        return 0
    end
end