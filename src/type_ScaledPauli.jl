"""
    coeff::T
    pauli::Pauli{N}
"""
struct ScaledPauli{T,N}
    coeff::T
    pauli::Pauli{N}
end

# ScaledPauliVector{T,N} = Vector{ScaledPauli{T,N}}

function Base.convert(::Type{ScaledPauli{T,N}}, p::Pauli{N}) where {T,N}
    return ScaledPauli{T,N}(get_phase(p), phase(phasefree(p)))
end

function ScaledPauli(p::Pauli{N}) where N
    return ScaledPauli{ComplexF64, N}(get_phase(p), phasefree(p)) 
end

function Base.display(sp::ScaledPauli)
    if sp.pauli.θ == 0
        @printf(" %12.8f %12.8fi   1  %s\n", real(sp.coeff), imag(sp.coeff), string(sp.pauli))
    elseif sp.pauli.θ == 1                                        
        @printf(" %12.8f %12.8fi   i  %s\n", real(sp.coeff), imag(sp.coeff), string(sp.pauli))
    elseif sp.pauli.θ == 2                                        
        @printf(" %12.8f %12.8fi  -1  %s\n", real(sp.coeff), imag(sp.coeff), string(sp.pauli))
    elseif sp.pauli.θ == 3                                        
        @printf(" %12.8f %12.8fi  -i  %s\n", real(sp.coeff), imag(sp.coeff), string(sp.pauli))
    end
end
function Base.display(sv::Vector{ScaledPauli{T,N}}) where {T,N}
    for i in sv
        display(i)
    end
end


"""
    Base.:-(p::ScaledPauli{T,N}) where {T,N}

Negate the `p`. We could either change the coeff or the pauli. Not sure which is best, negate coeff currently
"""
function Base.:-(p::ScaledPauli{T,N}) where {T,N}
    return ScaledPauli{T,N}(-p.coeff, p.pauli)
end
function Base.:-(p::Pauli{N}) where {N}
    return rotate_phase(p,2) 
end

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
    commute(p1::ScaledPauli{T,N}, p2::ScaledPauli{T,N}) where {T,N}

TBW
"""
function commute(p1::ScaledPauli{T,N}, p2::ScaledPauli{T,N}) where {T,N}
    return commute(p1.pauli, p2.pauli)
end
function commute(p1::Pauli{N}, p2::ScaledPauli{T,N}) where {T,N}
    return commute(p1, p2.pauli)
end
function commute(p1::ScaledPauli{T,N}, p2::Pauli{N}) where {T,N}
    return commute(p1.pauli, p2)
end


get_coeff(p::ScaledPauli{T,N}) where {T,N} = p.coeff * get_phase(p.pauli)
get_coeff(p::Pauli{N}) where {N} = get_phase(p)


Base.isless(p1::ScaledPauli{T,N}, p2::ScaledPauli{T,N}) where {T,N} = isless(p1.pauli, p2.pauli)
Base.isequal(p1::ScaledPauli{T,N}, p2::ScaledPauli{T,N}) where {T,N} = isequal(p1.pauli, p2.pauli)

function Base.unique!(spv::Vector{ScaledPauli{T,N}}) where {T,N}
    sort!(spv)

    # i = 1
    # while i < length(spv)
    #     njump = 0
    #     j = i+1
    #     while j < length(spv) && isequal(spv[j], spv[i])
    #         println(i, " ", j)
    #         spv[i] = ScaledPauli{T,N}(get_coeff(spv[i]) + get_coeff(spv[j]), phasefree(spv[i].pauli))
    #         spv[j] *= 0 
    #         j += 1
    #         njump += 1
    #     end
    #     i += 1 + njump
    #     j += 1
    #     # for j in i+1:length(spv)
    #     #     if isless(spv[i], spv[j])
    #     #         # spv[i+1] = spv[j]
    #     #         i += njump
    #     #         continue
    #     #     elseif isequal(spv[i], spv[j])
    #     #         njump += 1
    #     #         spv[i] = ScaledPauli{T,N}(get_coeff(spv[i]) + get_coeff(spv[j]), phasefree(spv[i].pauli))
    #     #     else
    #     #         error(i, " ", j, " Why is this not sorted?")
    #     #     end 
    #     # end
    # end
    for i in 1:length(spv)-1
        for j in i+1:length(spv)
            if isless(spv[i], spv[j])
                continue
            elseif isequal(spv[i], spv[j])
                
            end
        end
    end
end
