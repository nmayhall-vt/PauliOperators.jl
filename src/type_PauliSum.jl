"""
    ops::Dict{Pauli{N},ComplexF64}

A collection of `Pauli`s, joined by addition
"""
struct PauliSum{N}  
    ops::Dict{Pauli{N},ComplexF64}
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

Base.get(ps::PauliSum{N}, p::Pauli{N}) where N = get(ps.ops, p, zero(ComplexF64))
Base.keys(ps::PauliSum) = keys(ps.ops)
Base.getindex(ps::PauliSum{N}, p::Pauli{N}) where N = ps.ops[p]
Base.setindex!(ps::PauliSum{N}, v, p::Pauli{N}) where N = ps.ops[p] = v

"""
    Base.+(p1::PauliSum{N}, p2::PauliSum{N}) where {N}

Add two `PauliSum`s. 
"""
function Base.:+(ps1::PauliSum{N}, ps2::PauliSum{N}) where {N}
    return PauliSum{N}(mergewith(+, ps1.ops, ps2.ops))
end

"""
    Base.sum!(p1::PauliSum{N}, p2::PauliSum{N}) where {N}

Add two `PauliSum`s. 
"""
function Base.sum!(ps1::PauliSum{N}, ps2::PauliSum{N}) where {N}
    mergewith!(+, ps1.ops, ps2.ops)
end

"""
    Base.-(p1::PauliSum{N}, p2::PauliSum{N}) where {N}

Subtract two `PauliSum`s. 
"""
function Base.:-(ps1::PauliSum{N}, ps2::PauliSum{N}) where {N}
    return PauliSum{N}(mergewith(-, ps1.ops, ps2.ops))
end