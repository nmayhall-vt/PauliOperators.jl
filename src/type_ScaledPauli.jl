"""
    coeff::T
    pauli::Pauli{N}

Simply a combination of a `Pauli` with a coefficient. When sorted, only `pauli` is considered for the comparisons.
"""
struct ScaledPauli{N} <: AbstractPauli{N}
    coeff::ComplexF64
    pauli::FixedPhasePauli{N}
end



function ScaledPauli(p::Pauli{N}) where N
    return ScaledPauli{N}(get_phase(p), p.pauli) 
end

ScaledPauliVector{N} = Vector{ScaledPauli{N}}
function ScaledPauliVector(N)
    return Vector{ScaledPauli{N}}([])
end

function Base.display(sp::ScaledPauli)
    @printf(" %12.8f %12.8fi   %s\n", real(sp.coeff), imag(sp.coeff), string(sp.pauli))
end
function Base.display(sv::Vector{ScaledPauli{N}}) where {N}
    for i in sv
        display(i)
    end
end
"""
    Base.display(p::Pauli)

Display, y = iY
"""
Base.string(p::ScaledPauli) = string(p.pauli)


"""
    Base.:-(p::ScaledPauli{N}) where {T,N}

Negate the `p`. We could either change the coeff or the pauli. Not sure which is best, negate coeff currently
"""
function Base.:-(p::ScaledPauli{N}) where {N}
    return ScaledPauli{N}(-p.coeff, p.pauli)
end
Base.:≈(p::ScaledPauli{N}, q::ScaledPauli{N}) where {N} = (p.pauli == q.pauli) & (p.coeff ≈ q.coeff)




get_coeff(p::ScaledPauli{N}) where {T,N} = p.coeff
get_coeff(p::Pauli{N}) where {N} = get_phase(p)
# get_coeff(p::ScaledPauli{N}) where {T,N} = p.coeff * get_phase(p.pauli)
# get_coeff(p::Pauli{N}) where {N} = get_phase(p)


Base.isless(p1::ScaledPauli{N}, p2::ScaledPauli{N}) where {T,N} = isless(p1.pauli, p2.pauli)
Base.isequal(p1::ScaledPauli{N}, p2::ScaledPauli{N}) where {T,N} = isequal(p1.pauli, p2.pauli)

"""
    Base.unique!(spv::Vector{ScaledPauli{N}}) where {T,N}

Sort and add duplicates. This modifies `spv` to be a sorted list of only unique `ScaledPauli`'s
"""
function Base.unique!(spv::Vector{ScaledPauli{N}}) where {T,N}
    sort!(spv)
    i = 1
    for j in 2:length(spv)
        if isless(spv[i], spv[j])
            i += 1
            spv[i] = spv[j] 
        elseif isequal(spv[i], spv[j])
            spv[i] = ScaledPauli{N}(get_coeff(spv[i]) + get_coeff(spv[j]), phasefree(spv[i].pauli))
        end
    end
    resize!(spv,i)
end

"""
    Base.unique(spv::Vector{ScaledPauli{N}}) where {T,N}

Sort and add duplicates. This returns a sorted list of only unique `ScaledPauli`'s
"""
function Base.unique(spv::Vector{ScaledPauli{N}}) where {T,N}
    out = deepcopy(spv)
    unique!(out)
    return out
end


"""
    Base.Matrix(spv::Vector{ScaledPauli{N}}) where {T,N}

TBW
"""
function Base.Matrix(spv::Vector{ScaledPauli{N}}; T=ComplexF64) where {N}
    out = zeros(T, Int128(2)^N, Int128(2)^N)
    for spvi in spv 
        out .+= Matrix(spvi.pauli) .* spvi.coeff
    end
    return out
end


"""
    Base.Matrix(p::Pauli{N}) where N

Create dense matrix representation 
"""
function Base.Matrix(p::ScaledPauli{N}) where {T,N}
    return Matrix(p.pauli) .* p.coeff 
end


"""
    random_ScaledPauli(N)

TBW
"""
function random_ScaledPauli(N)
    return ScaledPauli{N}(rand(ComplexF64), random_FixedPhasePauli(N))
end

Base.adjoint(sp::ScaledPauli{N}) where {N} = is_hermitian(sp.pauli) ? ScaledPauli{N}(adjoint(sp.coeff), sp.pauli) : ScaledPauli{N}(-adjoint(sp.coeff), sp.pauli)