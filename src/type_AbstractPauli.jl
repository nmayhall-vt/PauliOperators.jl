abstract type AbstractPauli{N} end


@inline nqubits(p::AbstractPauli{N}) where N = N
@inline nY(p::AbstractPauli) = count_ones(p.x & p.z)




"""
    commute(p1::AbstractPauli, p2::AbstractPauli)

Check if they commute, return true or false
"""
function commute(p1::AbstractPauli, p2::AbstractPauli)
    return iseven(count_ones(p1.x & p2.z) - count_ones(p1.z & p2.x)) 
end