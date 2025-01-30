Base.:*(p::Pauli, pb::PauliBasis) = p*Pauli(pb)
Base.:*(pb::PauliBasis, p::Pauli) = Pauli(pb)*p
Base.:*(ps::PauliSum, p::Union{Pauli, PauliBasis}) = ps * PauliSum(p) 
Base.:*(p::Union{Pauli, PauliBasis}, ps::PauliSum) = PauliSum(p) * ps 

function Base.:*(d1::Union{Dyad{N}, DyadBasis{N}}, d2::Union{Dyad{N}, DyadBasis{N}}) where N
    return Dyad{N}(coeff(d1) * coeff(d2) * (d1.bra.v==d2.ket.v), d1.ket, d2.bra)
end

Base.:*(k::Ket{N}, b::Bra{N}) where N = DyadBasis{N}(k,b)
Base.:*(k::Bra{N}, b::Ket{N}) where N = k.v == b.v ? 1 : 0 

function Base.:*(b::Bra{N}, p::Pauli{N}) where N
    new_bra = Bra{N}(b.v ⊻ p.x)
    sign = count_ones(p.z & b.v)%2
    return sign == 0 ? p.s : -p.s, new_bra
end

function Base.:*(b::Bra{N}, p::PauliBasis{N}) where N
    new_bra = Bra{N}(b.v ⊻ p.x)
    sign = count_ones(b.v & p.z)%2
    return sign == 0 ? 1im^symplectic_phase(p) : -1im^symplectic_phase(p), new_bra
end

function Base.:*(p::Pauli{N}, k::Ket{N}) where N
    new_ket = Ket{N}(p.x ⊻ k.v)
    sign = count_ones(p.z & new_ket.v)%2
    return sign == 0 ? p.s : -p.s, new_ket
end

function Base.:*(p::PauliBasis{N}, k::Ket{N}) where N
    new_ket = Ket{N}(p.x ⊻ k.v)
    sign = count_ones(p.z & new_ket.v)%2
    return sign == 0 ? 1im^symplectic_phase(p) : -1im^symplectic_phase(p), new_ket
end

function Base.:*(p::Union{Pauli{N}, PauliBasis{N}}, d::Union{Dyad{N}, DyadBasis{N}}) where N
    new_coeff, new_ket = p*d.ket
    return Dyad{N}(new_coeff * coeff(d) , new_ket, d.bra)
end 
function Base.:*(d::Union{Dyad{N}, DyadBasis{N}}, p::Union{Pauli{N}, PauliBasis{N}}) where N
    new_coeff, new_bra = d.bra*p
    return Dyad{N}(new_coeff * coeff(d) , d.ket, new_bra)
end 