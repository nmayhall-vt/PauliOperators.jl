Base.isless(p1::Pauli{N}, p2::Pauli{N}) where N = isless((p1.z, p1.x), (p2.z, p2.x))
Base.isless(p1::ScaledPauli{T,N}, p2::Pauli{N}) where {T,N} = isless(p1.pauli, p2)
Base.isless(p1::Pauli{N}, p2::ScaledPauli{T,N}) where {T,N} = isless(p1, p2.pauli)

Base.adjoint(p::Pauli{N}) where N = is_hermitian(p) ? p : -p
Base.adjoint(sp::ScaledPauli{T,N}) where {T,N} = is_hermitian(sp.pauli) ? ScaledPauli{T,N}(adjoint(sp.coeff), sp.pauli) : ScaledPauli{T,N}(-adjoint(sp.coeff), sp.pauli)


# function Base.adjoint(spv::Vector{ScaledPauli{T,N}}) where {T,N}
#     # for i in 1:length(spv)
#     #     spv[i] = adjoint(spv[i])
#     # end
# end






####################################################
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

"""
    Base.:+(ps::PauliSum{N}, p::Pauli{N}) where {N}

Add a `Pauli` to a PauliSum. 
"""
function Base.:+(ps::PauliSum{N}, p::Pauli{N}) where {N}
    out = deepcopy(ps)
    sum!(out, p)
    return out
end
Base.:+(p::Pauli{N}, ps::PauliSum{N}) where {N} = ps + p


"""
    Base.sum!(p1::PauliSum{N}, p2::Pauli{N}) where {N}

Add a `Pauli` to a PauliSum. 
"""
Base.sum!(ps::PauliSum{N}, p::Pauli{N}) where {N} = ps[phasefree(p)] = get(ps, p) + get_phase(p) 
Base.sum!(ps::Vector{ScaledPauli{T,N}}, p::ScaledPauli{T,N}) where {T,N} = push!(ps, p)
Base.sum!(ps::Vector{ScaledPauli{T,N}}, p::Pauli{N}) where {T,N} = push!(ps, ScaledPauli(p))

function Base.:+(p1::Vector{ScaledPauli{T,N}}, p2::Vector{ScaledPauli{T,N}}) where {T,N}
    out = deepcopy(p1)
    append!(out, p2)
    return out 
end
Base.:+(p::Pauli{N}, ps::PauliSum{N}) where {N} = ps + p
