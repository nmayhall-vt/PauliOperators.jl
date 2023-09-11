abstract type AbstractPauli{N} end


@inline nqubits(p::AbstractPauli{N}) where N = N
@inline nY(p::AbstractPauli) = count_ones(p.x & p.z)



