# 
function Base.:+(p1::ScaledPauli{N}, p2::ScaledPauli{N}) where {N}
    if isequal(p1.pauli, p2.pauli)
        return PauliSum{N}(Dict{FixedPhasePauli{N},ComplexF64}(p1.pauli=>p1.coeff+p2.coeff))
    else
        return PauliSum{N}(Dict{FixedPhasePauli{N},ComplexF64}(p1.pauli=>p1.coeff, p2.pauli=>p2.coeff))
    end
end
Base.:+(p::ScaledPauli{N}, a::Union{Pauli{N}, FixedPhasePauli{N}}) where {N} = p + ScaledPauli(a)
Base.:+(a::Union{Pauli{N}, FixedPhasePauli{N}}, p::ScaledPauli{N}) where {N} = ScaledPauli(a) + p
Base.:+(a::Union{Pauli{N}, FixedPhasePauli{N}}, b::Union{Pauli{N}, FixedPhasePauli{N}}) where {N} = ScaledPauli(a) + ScaledPauli(b)


Base.:+(p::ScaledPauli{N}, a::Number) where {N} = p + ScaledPauli{N}(a, FixedPhasePauli{N}(0,0))

function Base.:+(ps::PauliSum{N}, p::AbstractPauli{N}) where {N}
    out = deepcopy(ps)
    sum!(out, p)
    return out
end
Base.:+(p::AbstractPauli{N}, ps::PauliSum{N}) where {N} = ps + p

Base.:+(ps1::PauliSum{N}, ps2::PauliSum{N}) where {N} = PauliSum{N}(mergewith(+, ps1.ops, ps2.ops))


"""
    Base.sum!(p1::PauliSum{N}, p2::Pauli{N}) where {N}

Add a `Pauli` to a PauliSum. 
"""
Base.sum!(ps::PauliSum{N}, p::FixedPhasePauli{N}) where {N} = ps[p] = get(ps, p) + 1 
Base.sum!(ps::PauliSum{N}, p::Pauli{N}) where {N} = ps[p.pauli] = get(ps, p.pauli) + get_phase(p) 
Base.sum!(ps::PauliSum{N}, p::ScaledPauli{N}) where {N} = ps[p.pauli] = get(ps, p.pauli) + p.coeff 
Base.sum!(ps::Vector{ScaledPauli{N}}, p::ScaledPauli{N}) where {N} = push!(ps, p)
Base.sum!(ps::Vector{ScaledPauli{N}}, p::Pauli{N}) where {N} = push!(ps, ScaledPauli(p))
Base.sum!(ps::Vector{ScaledPauli{N}}, p::FixedPhasePauli{N}) where {N} = push!(ps, ScaledPauli(p))

function Base.:+(p1::Vector{ScaledPauli{N}}, p2::Vector{ScaledPauli{N}}) where {N}
    out = deepcopy(p1)
    append!(out, p2)
    return out 
end







Base.sum!(vs::SparseKetBasis{N,T}, k::KetBitString{N}) where {N,T} = vs[k] = get(vs, k, 0)+1 
Base.sum!(vs::SparseKetBasis{N,T}, k::KetBitString{N}, c) where {N,T} = vs[k] = get(vs, k, 0)+c 


