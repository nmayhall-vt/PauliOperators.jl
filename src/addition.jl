"""
    Base.sum!(p1::PauliSum{N}, p2::PauliSum{N}) where {N}

Add two `PauliSum`s. 
"""
function Base.sum!(ps1::PauliSum{N}, ps2::PauliSum{N}) where {N}
    mergewith!(+, ps1, ps2)
end

function Base.sum!(ps1::PauliSum{N,T}, ps2::Adjoint{<:Any, PauliSum{N,T}}) where {N,T}
    for (pauli, coeff) in ps2.parent
        if haskey(ps1, pauli)
            ps1[pauli] += coeff'
        else
            ps1[pauli] = coeff'
        end
    end
    return ps1
end
Base.:+(ps2::Adjoint{<:Any, PauliSum{N,T}}, ps1::PauliSum{N,T}) where {N,T} = ps1 + ps2 


function Base.:+(ps1::PauliSum, ps2::PauliSum)
    out = deepcopy(ps2)
    mergewith!(+, out, ps1)
    return out 
end

function Base.:+(ps1::DyadSum, ps2::DyadSum)
    out = deepcopy(ps2)
    mergewith!(+, out, ps1)
    return out 
end

"""
    Base.:-(ps1::PauliSum, ps2::PauliSum)

Subtract two `PauliSum`s. 
"""
function Base.:-(ps1::PauliSum, ps2::PauliSum)
    out = deepcopy(ps2)
    map!(x->-x, values(out))
    mergewith!(+, out, ps1)
    return out 
end


function Base.:-(ps1::DyadSum, ps2::DyadSum)
    out = deepcopy(ps2)
    map!(x->-x, values(out))
    mergewith!(+, out, ps1)
    return out 
end

function Base.:+(ps1::PauliSum{N,T}, ps2::Adjoint{<:Any, PauliSum{N,T}}) where {N,T}
    out = deepcopy(ps1)
    sum!(out, ps2)
    return out
end

Base.sum!(a::PauliSum{N}, b::Union{Pauli{N}, PauliBasis{N}}) where N = sum!(a, promote_to_sum(b)) 
Base.sum!(a::DyadSum{N}, b::Union{Dyad{N}, DyadBasis{N}}) where N = sum!(a, promote_to_sum(b)) 
Base.sum!(a::Union{Pauli{N}, PauliBasis{N}}, b::PauliSum{N}) where N = sum!(b, promote_to_sum(a)) 
Base.sum!(a::Union{Dyad{N}, DyadBasis{N}}, b::DyadSum{N}) where N = sum!(b, promote_to_sum(a)) 


Base.:+(a::Singles{N}, b::Singles{N}) where N = promote_to_sum(a) + promote_to_sum(b) 
Base.:-(a::Singles{N}, b::Singles{N}) where N = promote_to_sum(a) - promote_to_sum(b) 
Base.:+(a::Singles{N}, b::Sums{N,T}) where {N,T} = promote_to_sum(a) + b
Base.:+(a::Sums{N,T}, b::Singles{N}) where {N,T} = a + promote_to_sum(b)
Base.:-(a::Singles{N}, b::Sums{N,T}) where {N,T} = promote_to_sum(a) - b
Base.:-(a::Sums{N,T}, b::Singles{N}) where {N,T} = a - promote_to_sum(b)

