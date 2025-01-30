Base.:*(p::Pauli, pb::PauliBasis) = p*Pauli(pb)
Base.:*(pb::PauliBasis, p::Pauli) = Pauli(pb)*p
Base.:*(ps::PauliSum, p::Union{Pauli, PauliBasis}) = ps * PauliSum(p) 
Base.:*(p::Union{Pauli, PauliBasis}, ps::PauliSum) = PauliSum(p) * ps 

Base.:*(k::Ket{N}, b::Bra{N}) where N = DyadBasis{N}(k,b)
Base.:*(k::Bra{N}, b::Ket{N}) where N = k.v == b.v ? 1 : 0 

