
"""
    index(k::KetBitString{N})

Return location of this basis vector in the basis of `KetBitStrings`
"""
index(k::KetBitString) = 1 + k.v

"""
    index(d::Dyad{N}) where N

Return location of this basis vector in the basis of `Dyad`'s
"""
index(d::Dyad{N}) where N = 1 + d.ket.v + d.bra.v*(BigInt(2)^N)

"""
    index(p::FixedPhasePauli{N}) where N

Return location of this basis vector in the basis of `FixedPhasePauli`'s
"""
index(p::FixedPhasePauli{N}) where N = 1 + p.z + p.x*(BigInt(2)^N)