abstract type AbstractPauli{N} end

# abstract type BasisTrait end
# struct IsBasisVector <: BasisTrait end
# struct NotBasisVector <: BasisTrait end

# BasisTrait(::Type) = NotBasisVector()
# BasisTrait(::FixedPhasePauli{N}) where N = IsBasisVector()
# BasisTrait(::Dyad{N}) where N = IsBasisVector()
# BasisTrait(::KetBitString{N}) where N = IsBasisVector()


@inline nqubits(p::AbstractPauli{N}) where N = N
@inline nY(p::AbstractPauli) = count_ones(p.x & p.z)



