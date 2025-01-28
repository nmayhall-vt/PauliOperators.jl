Base.:*(p::Pauli, pb::PauliBasis) = p*Pauli(pb)
Base.:*(pb::PauliBasis, p::Pauli) = Pauli(pb)*p
Base.:*(ps::PauliSum, p::Union{Pauli, PauliBasis}) = ps * PauliSum(p) 
Base.:*(p::Union{Pauli, PauliBasis}, ps::PauliSum) = PauliSum(p) * ps 