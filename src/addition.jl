function Base.:+(ps::PauliSum{N, T}, p::Union{Pauli{N}, PauliBasis{N}}) where {N, T}
    out = deepcopy(ps)
    sum!(out, p)
    return out
end

function Base.sum!(ps::PauliSum{N, T}, p::Union{Pauli{N}, PauliBasis{N}}) where {N,T}
    if haskey(ps, PauliBasis(p))
        ps[PauliBasis(p)] += coeff(p) 
    else
        ps[PauliBasis(p)] = coeff(p)
    end
    return ps
end
Base.sum!(p::Union{Pauli{N}, PauliBasis{N}}, ps::PauliSum{N, T}) where {N,T} = sum!(ps, p)
Base.:+(p::Union{Pauli{N}, PauliBasis{N}}, ps::PauliSum{N, T}) where {N,T} = ps + p



function Base.sum!(ps::DyadSum{N, T}, p::Union{Dyad{N}, DyadBasis{N}}) where {N,T}
    if haskey(ps, DyadBasis(p))
        ps[DyadBasis(p)] += coeff(p) 
    else
        ps[DyadBasis(p)] = coeff(p)
    end
    return ps
end

function Base.:+(ps::DyadSum{N, T}, p::Union{Dyad{N}, DyadBasis{N}}) where {N, T}
    out = deepcopy(ps)
    sum!(out, p)
    return out
end

Base.sum!(p::Union{Dyad{N}, DyadBasis{N}}, ps::DyadSum{N, T}) where {N,T} = sum!(ps,p)
Base.:+(p::Union{Dyad{N}, DyadBasis{N}}, ps::DyadSum{N, T}) where {N,T} = ps + p