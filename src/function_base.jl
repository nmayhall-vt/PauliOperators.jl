Base.isless(p1::FixedPhasePauli{N}, p2::FixedPhasePauli{N}) where N = isless((p1.z, p1.x), (p2.z, p2.x))
Base.isless(p1::Pauli{N}, p2::Pauli{N}) where N = isless((p1.z, p1.x), (p2.z, p2.x))
Base.isless(p1::ScaledPauli{N}, p2::Pauli{N}) where {N} = isless(p1.pauli, p2)
Base.isless(p1::Pauli{N}, p2::ScaledPauli{N}) where {N} = isless(p1, p2.pauli)

# Base.adjoint(p::Pauli{N}) where N = is_hermitian(p) ? p : -p
# Base.adjoint(sp::ScaledPauli{N}) where {T,N} = is_hermitian(sp.pauli) ? ScaledPauli{N}(adjoint(sp.coeff), sp.pauli) : ScaledPauli{N}(-adjoint(sp.coeff), sp.pauli)


# # function Base.adjoint(spv::Vector{ScaledPauli{N}}) where {T,N}
# #     # for i in 1:length(spv)
# #     #     spv[i] = adjoint(spv[i])
# #     # end
# # end






# ####################################################
