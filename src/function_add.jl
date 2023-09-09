# 
function Base.:+(p1::ScaledPauli{T,N}, p2::ScaledPauli{T,N}) where {T,N}
    if isequal(p1.pauli, p2.pauli)
        return Vector{ScaledPauli{T,N}}([ScaledPauli{T,N}(get_coeff(p1) + get_coeff(p2), phasefree(p1.pauli))])
    else
        return Vector{ScaledPauli{T,N}}([p1, p2])
    end
end
function Base.:+(p::ScaledPauli{T,N}, a::Number) where {T,N}
    return p + ScaledPauli{T,N}(a, Pauli{N}(0,0,0))
end
function Base.:+(p::ScaledPauli{T,N}, a::Pauli{N}) where {T,N}
    return p + ScaledPauli{T,N}(1, a)
end

"""
    Base.:+(p1::Pauli{N}, p2::Pauli{N}) where {N}

Add two `Pauli`'s together. This returns a `PauliSum`
"""
function Base.:+(p1::Pauli{N}, p2::Pauli{N}) where {N}
    if phasefree(p1) == phasefree(p2) 
        return PauliSum{N}(Dict(phasefree(p1)=>get_phase(p1)+get_phase(p2)))
    else
        return PauliSum{N}(Dict(phasefree(p1)=>get_phase(p1), phasefree(p2)=>get_phase(p2)))
    end
end
function Base.:+(p1::Pauli{N}, p2::PauliPF{N}) where {N}
    if p1.z == p2.z & p1.x == p2.x 
        return PauliSum{N}(Dict(phasefree(p1)=>get_phase(p1)+get_phase(p2)))
    else
        return PauliSum{N}(Dict(phasefree(p1)=>get_phase(p1), Pauli(p2)=>get_phase(p2)))
    end
end
function Base.:+(p1::PauliPF{N}, p2::PauliPF{N}) where {N}
    if p1 == p2 
        return PauliSum{N}(Dict(p1=>2))
    else
        return PauliSum{N}(Dict(p1=>1, p2=>1))
    end
end
function Base.:+(p1::Pauli{N}, p2::ScaledPauli{T,N}) where {T,N}
    if phasefree(p1) == phasefree(p2) 
        return PauliSum{N}(Dict(phasefree(p1)=>get_phase(p1)+get_phase(p2.pauli)*p2.coeff))
    else
        return PauliSum{N}(Dict(phasefree(p1)=>get_phase(p1), phasefree(p2)=>get_phase(p2.pauli)*p2.coeff))
    end
end
function Base.:+(p1::PauliPF{N}, p2::ScaledPauli{T,N}) where {T,N}
    if p1 == p2.pauli 
        return PauliSum{N}(Dict(p1 => 1 + p2.coeff))
    else
        return PauliSum{N}(Dict(p1 => 1, p2.pauli => p2.coeff))
    end
end

Base.:+(p1::ScaledPauli{T,N}, p2::Pauli{N}) where {T,N} = p2 + p1
Base.:+(p1::ScaledPauli{T,N}, p2::PauliPF{N}) where {T,N} = p2 + p1

function Base.:+(ps::PauliSum{N}, p::Pauli{N}) where {N}
    out = deepcopy(ps)
    sum!(out, p)
    return out
end
function Base.:+(ps::PauliSum{N}, p::PauliPF{N}) where {N}
    out = deepcopy(ps)
    sum!(out, p)
    return out
end
Base.:+(p::Pauli{N}, ps::PauliSum{N}) where {N} = ps + p
Base.:+(p::PauliPF{N}, ps::Pauli{N}) where {N} = ps + p
Base.:+(ps1::PauliSum{N}, ps2::PauliSum{N}) where {N} = PauliSum{N}(mergewith(+, ps1.ops, ps2.ops))


"""
    Base.sum!(p1::PauliSum{N}, p2::Pauli{N}) where {N}

Add a `Pauli` to a PauliSum. 
"""
Base.sum!(ps::PauliSum{N}, p::Pauli{N}) where {N} = ps[phasefree(p)] = get(ps, p) + get_phase(p) 
Base.sum!(ps::PauliSum{N}, p::PauliPF{N}) where {N} = ps[p] = get(ps, p) + 1 
Base.sum!(ps::Vector{ScaledPauli{T,N}}, p::ScaledPauli{T,N}) where {T,N} = push!(ps, p)
Base.sum!(ps::Vector{ScaledPauli{T,N}}, p::Pauli{N}) where {T,N} = push!(ps, ScaledPauli(p))
Base.sum!(ps::Vector{ScaledPauli{T,N}}, p::PauliPF{N}) where {T,N} = push!(ps, ScaledPauli(p))

function Base.:+(p1::Vector{ScaledPauli{T,N}}, p2::Vector{ScaledPauli{T,N}}) where {T,N}
    out = deepcopy(p1)
    append!(out, p2)
    return out 
end






