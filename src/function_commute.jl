"""
    commute(p1,p2)

Check if `p1` and `p2` commute, return true or false
"""
@inline commute(p1::FixedPhasePauli, p2::FixedPhasePauli) = iseven(count_ones(p1.x & p2.z) - count_ones(p1.z & p2.x)) 
@inline commute(p1::Pauli{N}, p2::Pauli{N}) where {N} = commute(p1.pauli, p2.pauli)
@inline commute(p1::ScaledPauli{N}, p2::ScaledPauli{N}) where {N} = commute(p1.pauli, p2.pauli)
@inline commute(p1::Pauli{N}, p2::ScaledPauli{N}) where {N} = commute(p1, p2.pauli)
@inline commute(p1::ScaledPauli{N}, p2::Pauli{N}) where {N} = commute(p1.pauli, p2)

function commutator(spA::ScaledPauli{N}, spB::ScaledPauli{N}) where {N}
    res = ScaledPauli{N}[]
    if commute(spA.pauli, spB.pauli)
        return res
    else
        comm = ScaledPauli{N}[]
        comm1 = spA*spB; comm2 = -spB*spA;
        push!(comm, comm1); push!(comm, comm2)
        res = unique(comm)
    end
    unique!(res); filter!(sp->abs(sp.coeff) > 1.0e-10, res)
    return res    # return type: Vector{ScaledPauli{N}}
end

function commutator(spv1::Vector{ScaledPauli{N}}, spv2::Vector{ScaledPauli{N}}) where {N}
    res = ScaledPauli{N}[]
    for sp1 in spv1
        for sp2 in spv2
            comm = commutator(sp1,sp2)
            if !isempty(comm)
                for sp in comm
                    push!(res, sp)
                end
            end
        end
    end
    unique!(res); filter!(sp->abs(sp.coeff) > 1.0e-10, res)
    return res      # return type: Vector{Vector{ScaledPauli{N}}}
end


# function commutator(sdA::ScaledDyad{N,T}, sdB::ScaledDyad{N,T}) where {N,T}
#     res = ScaledDyad{N}[]
#     if commute(spA.pauli, spB.pauli)
#         return res
#     else
#         comm = ScaledPauli{N}[]
#         comm1 = spA*spB; comm2 = -spB*spA;
#         push!(comm, comm1); push!(comm, comm2)
#         res = unique(comm)
#     end
#     unique!(res); filter!(sp->abs(sp.coeff) > 1.0e-10, res)
#     return res    # return type: Vector{ScaledPauli{N}}
# end
