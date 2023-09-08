
"""
    otimes(p1::Pauli{N}, p2::Pauli{M}) where {N,M}

TBW
"""
function otimes(p1::Pauli{N}, p2::Pauli{M}) where {N,M}
    return Pauli{N+M}((p1.θ + p2.θ)%4, p1.z | p2.z << N, p1.x | p2.x << N)
end
const ⊗ = otimes

"""
    otimes(p1::PauliSum{N}, p2::Pauli{M}) where {N,M}

TBW
"""
function otimes(p1::PauliSum{N}, p2::Pauli{M}) where {N,M}
    out = PauliSum(N+M)
    for (op,coeff) in p1.ops
        out.ops[phasefree(op ⊗ p2)] = coeff * get_phase(p2)
    end
    return out 
end
"""
    otimes(p::Pauli{N}, ps::PauliSum{M}) where {N,M}

TBW
"""
function otimes(p::Pauli{N}, ps::PauliSum{M}) where {N,M}
    out = PauliSum(N+M)
    for (op,coeff) in ps.ops
        out.ops[phasefree(p ⊗ op)] = coeff * get_phase(p)
    end
    return out 
end

"""
    otimes(p1::PauliSum{N}, p2::PauliSum{M}) where {N,M}

TBW
"""
function otimes(p1::PauliSum{N}, p2::PauliSum{M}) where {N,M}
    out = PauliSum(N+M)
    for (op1,coeff1) in p1.ops
        for (op2,coeff2) in p2.ops
            out.ops[phasefree(op1 ⊗ op2)] = coeff1 * coeff2 * get_phase(op1) * get_phase(op2)
        end
    end
    return out 
end

"""
    osum(p1::Pauli{N}, p2::Pauli{M}) where {N,M}

Returns the direct sum of two Paulis
"""
function osum(p1::Pauli{N}, p2::Pauli{M}) where {N,M}
    return p1 ⊗ Pauli(M) + Pauli(N) ⊗ p2 
end
const ⊕ = osum

"""
    osum(p1::PauliSum{N}, p2::Pauli{M}) where {N,M}

Returns the direct sum of a PauliSum and a Pauli
"""
function osum(p1::PauliSum{N}, p2::Pauli{M}) where {N,M}
    return p1 ⊗ Pauli(M) + Pauli(N) ⊗ p2 
end

"""
    osum(p1::Pauli{N}, p2::PauliSum{M}) where {N,M}

Returns the direct sum of a PauliSum and a Pauli
"""
function osum(p1::Pauli{N}, p2::PauliSum{M}) where {N,M}
    return p1 ⊗ Pauli(M) + Pauli(N) ⊗ p2 
end

"""
    osum(p1::Pauli{N}, p2::Pauli{M}) where {N,M}

Returns the direct sum of two PauliSums
"""
function osum(p1::PauliSum{N}, p2::PauliSum{M}) where {N,M}
    return p1 ⊗ Pauli(M) + Pauli(N) ⊗ p2 
end

"""
    Base.:(==)(p1::Pauli{N}, p2::Pauli{N}) where {N}

Check if they are equal, return true or false
"""
function Base.:(==)(p1::Pauli{N}, p2::Pauli{N}) where {N}
    return p1.x == p2.x && p1.z == p2.z && p1.θ == p2.θ
end
function Base.isequal(p1::Pauli{N}, p2::Pauli{N}) where N
    return p1.x == p2.x && p1.z == p2.z
end


"""
    Base.:(>)(p1::Pauli{N}, p2::Pauli{N}) where {N}

Check if `p1` > `p2`
"""
function Base.:>(p1::Pauli{N}, p2::Pauli{N}) where {N}
    return p1.z > p2.z || p1.z == p2.z && p1.x > p2.x
end


"""
    Base.:(<)(p1::Pauli{N}, p2::Pauli{N}) where {N}

Check if `p1` < `p2`
"""
function Base.:<(p1::Pauli{N}, p2::Pauli{N}) where {N}
    return p1.z < p2.z || p1.z == p2.z && p1.x < p2.x
end




# """
#     Base.sum!(p1::PauliSum{N}, p2::Pauli{N}) where {N}

# Add a `Pauli` to a PauliSum. 
# """
# function subtract!(p::Pauli{N}, ps::PauliSum{N}) where {N}
#     ps[phasefree(p)] = get_phase(p) - get(ps, p)
# end

# """
#     Base.sum!(p1::PauliSum{N}, p2::Pauli{N}) where {N}

# Add a `Pauli` to a PauliSum. 
# """
# function subtract!(ps::PauliSum{N}, p::Pauli{N}) where {N}
#     ps[phasefree(p)] = get(ps, p) - get_phase(p) 
# end

# """
#     Base.:-(p1::Pauli{N}, p2::Pauli{N}) where {N}

# Subtract two `Pauli`'s. This returns a `PauliSum`
# """
# function Base.:-(p1::Pauli{N}, p2::Pauli{N}) where {N}
#     if Base.isequal(p1, p2)
#         return PauliSum{N}(Dict(phasefree(p1)=>get_phase(p1)-get_phase(p2)))
#     else
#         return PauliSum{N}(Dict(phasefree(p1)=>get_phase(p1), phasefree(p2)=>-get_phase(p2)))
#     end
# end

# """
#     Base.:-(ps::PauliSum{N}, p::Pauli{N}) where {N}

# Add a `Pauli` to a PauliSum. 
# """
# function Base.:-(ps::PauliSum{N}, p::Pauli{N}) where {N}
#     out = deepcopy(ps)
#     subtract!(out, p)
#     return out
# end

# """
#     Base.:-(ps::PauliSum{N}, p::Pauli{N}) where {N}

# Add a `Pauli` to a PauliSum. 
# """
# function Base.:-(p::Pauli{N}, ps::PauliSum{N}) where {N}
#     out = deepcopy(ps)
#     subtract!(p, out)
#     return out
# end

"""
    Base.hash(p::Pauli{N}, h::UInt) where N

Create a hash for a `Pauli`. Because we want to collect matching operators, 
    with different phases, we don't actually put the phase in the hash
"""
function Base.hash(p::Pauli{N}, h::UInt) where N
    return hash((p.θ, p.z, p.x), h)
end
function Base.hash(p::Pauli{N}) where N
    return hash((p.θ, p.z, p.x))
end

"""
    get_phase(p::Pauli)

Return the phase of the `Pauli`, i^θ
"""
function get_phase(p::Pauli)
    return 1im^p.θ
end

"""
    rotate_phase(p::Pauli{N}, θ::Integer) where N

Rotate phase in units of π/2. In otherwords, multiply the phase by i^θ.
E.g., mutliplication by -1 is obtained with θ=2.
"""
function rotate_phase(p::Pauli{N}, θ::Integer) where N
    return Pauli{N}((p.θ + θ)%4, p.z, p.x)
end


"""
    negate(p::Pauli)

Multiply `p` by -1
"""
function negate(p::Pauli{N}) where N
    return rotate_phase(p,2) 
end

"""
    is_diagonal(p::Pauli)

Check if operator is diagonal in the computational (z) basis. E.g., does this operator consist of only I and/or Z?
"""
function is_diagonal(p::Pauli)
    return count_ones(p.x) == 0
end








"""
    expectation_value_sign(p::Pauli{N}, ket::Vector{Bool}) where N

compute expectation value of Pauli `o` for a product state `ket`
"""
function expectation_value_sign(p::Pauli{N}, ket::KetBitString{N}) where N
    is_diagonal(p) || return 0.0
    
    println(p)

    count_ones(p.z & ket.v) % 2 == 0 || return -(1im)^p.θ
    return (1im)^p.θ 
end
