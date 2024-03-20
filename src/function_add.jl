# 
function Base.:+(p1::ScaledPauli{N}, p2::ScaledPauli{N}) where {T,N}
    # if isequal(p1.pauli, p2.pauli)
    #     return Vector{ScaledPauli{N}}([ScaledPauli{N}(get_coeff(p1) + get_coeff(p2), phasefree(p1.pauli))])
    # else
    #     return Vector{ScaledPauli{N}}([p1, p2])
    # end
    if isequal(p1.pauli, p2.pauli)
        return PauliSum{N}(OrderedDict{FixedPhasePauli{N},ComplexF64}(p1.pauli=>p1.coeff+p2.coeff))
    else
        return PauliSum{N}(OrderedDict{FixedPhasePauli{N},ComplexF64}(p1.pauli=>p1.coeff, p2.pauli=>p2.coeff))
    end
end
function Base.:+(p::ScaledPauli{N}, a::Number) where {T,N}
    return p + ScaledPauli{N}(a, Pauli{N}(0,FixedPhasePauli{N}(0,0)))
end
Base.:+(p::ScaledPauli{N}, a::Pauli{N}) where {T,N} = p + ScaledPauli{N}(1, a.pauli)
Base.:+(a::Pauli{N}, p::ScaledPauli{N}) where {T,N} = p + a 

"""
    Base.:+(p1::Pauli{N}, p2::Pauli{N}) where {N}

Add two `Pauli`'s together. This returns a `PauliSum`
"""
function Base.:+(p1::Pauli{N}, p2::Pauli{N}) where {N}
    if phasefree(p1) == phasefree(p2) 
        return PauliSum{N}(OrderedDict(phasefree(p1)=>get_phase(p1)+get_phase(p2)))
    else
        return PauliSum{N}(OrderedDict(phasefree(p1)=>get_phase(p1), phasefree(p2)=>get_phase(p2)))
    end
end
function Base.:+(p1::Pauli{N}, p2::FixedPhasePauli{N}) where {N}
    if p1.z == p2.z & p1.x == p2.x 
        return PauliSum{N}(OrderedDict(phasefree(p1)=>get_phase(p1)+get_phase(p2)))
    else
        return PauliSum{N}(OrderedDict(phasefree(p1)=>get_phase(p1), Pauli(p2)=>get_phase(p2)))
    end
end
function Base.:+(p1::FixedPhasePauli{N}, p2::FixedPhasePauli{N}) where {N}
    if p1 == p2 
        return PauliSum{N}(OrderedDict(p1=>2))
    else
        return PauliSum{N}(OrderedDict(p1=>1, p2=>1))
    end
end


function Base.:+(ps::PauliSum{N}, p::FixedPhasePauli{N}) where {N}
    out = deepcopy(ps)
    sum!(out, p)
    return out
end
function Base.:+(ps::PauliSum{N}, p::Pauli{N}) where {N}
    out = deepcopy(ps)
    sum!(out, p)
    return out
end
function Base.:+(ps::PauliSum{N}, p::ScaledPauli{N}) where {N}
    out = deepcopy(ps)
    sum!(out, p)
    return out
end
Base.:+(p::FixedPhasePauli{N}, ps::PauliSum{N}) where {N} = ps + p
Base.:+(p::Pauli{N}, ps::PauliSum{N}) where {N} = ps + p
Base.:+(p::ScaledPauli{N}, ps::PauliSum{N}) where {N} = ps + p
Base.:+(p::FixedPhasePauli{N}, ps::Pauli{N}) where {N} = ps + p
# Base.:+(ps1::PauliSum{N}, ps2::PauliSum{N}) where {N} = PauliSum{N}(mergewith(+, ps1.ops, ps2.ops))
function Base.:+(ps1::PauliSum{N}, ps2::PauliSum{N}) where {N}
    out = deepcopy(ps1)
    mergewith!(+, out.ops, ps2.ops)
    return out 
end


"""
    Base.sum!(p1::PauliSum{N}, p2::Pauli{N}) where {N}

Add a `Pauli` to a PauliSum. 
"""
Base.sum!(ps::PauliSum{N}, p::FixedPhasePauli{N}) where {N} = ps[p] = get(ps, p) + 1 
Base.sum!(ps::PauliSum{N}, p::Pauli{N}) where {N} = ps[p.pauli] = get(ps, p.pauli) + get_phase(p) 
Base.sum!(ps::PauliSum{N}, p::ScaledPauli{N}) where {N} = ps[p.pauli] = get(ps, p.pauli) + p.coeff 
Base.sum!(ps::Vector{ScaledPauli{N}}, p::ScaledPauli{N}) where {T,N} = push!(ps, p)
Base.sum!(ps::Vector{ScaledPauli{N}}, p::Pauli{N}) where {T,N} = push!(ps, ScaledPauli(p))
Base.sum!(ps::Vector{ScaledPauli{N}}, p::FixedPhasePauli{N}) where {T,N} = push!(ps, ScaledPauli(p))
Base.sum!(ps::PauliSum{N}, p::PauliSum{N}) where {N} = PauliSum{N}(mergewith!(+, ps.ops, p.ops))

function Base.:+(p1::Vector{ScaledPauli{N}}, p2::Vector{ScaledPauli{N}}) where {T,N}
    out = deepcopy(p1)
    append!(out, p2)
    return out 
end







Base.sum!(vs::SparseKetBasis{N,T}, k::KetBitString{N}) where {N,T} = vs[k] = get(vs, k, 0)+1 
Base.sum!(vs::SparseKetBasis{N,T}, k::KetBitString{N}, c) where {N,T} = vs[k] = get(vs, k, 0)+c 