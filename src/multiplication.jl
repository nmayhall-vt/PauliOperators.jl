Base.:*(p::Pauli, pb::PauliBasis) = p*Pauli(pb)
Base.:*(pb::PauliBasis, p::Pauli) = Pauli(pb)*p
# Base.:*(ps::PauliSum, p::Union{Pauli, PauliBasis}) = ps * PauliSum(p) 
# Base.:*(p::Union{Pauli, PauliBasis}, ps::PauliSum) = PauliSum(p) * ps 
# Base.:*(p::Union{Pauli, PauliBasis}, ps::Adjoint{<:Any, PauliSum{N,T}}) where {N,T} = PauliSum(p) * ps 
# Base.:*(ps::Adjoint{<:Any, PauliSum{N,T}}, p::Union{Pauli, PauliBasis}) where {N,T} = ps * PauliSum(p)  

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


# function Base.:*(p::Union{Pauli{N}, PauliBasis{N}}, d::DyadSum{N,T}) where {N,T}
#     out = DyadSum(N)
#     for (dyad, coeff) in d
#         new_dyad = p*dyad
#         sum!(out, new_dyad * coeff)
#     end
#     return out 
# end 


# function Base.:*(d::DyadSum{N,T}, p::Union{Pauli{N}, PauliBasis{N}}) where {N,T}
#     out = DyadSum(N,T)
#     for (dyad, coeff) in d
#         new_dyad = dyad*p
#         sum!(out, new_dyad * coeff)
#     end
#     return out 
# end 

function Base.:*(d::DyadSum{N,T}, p::PauliSum{N,T}) where {N,T}
    out = DyadSum(N,T)
    for (dyad, coeff_d) in d
        for (pauli, coeff_p) in p
            new_dyad = dyad*pauli
            sum!(out, new_dyad * coeff_d * coeff_p)
        end   
    end
    return out 
end 

function Base.:*(d::DyadSum{N,T}, p::Adjoint{<:Any, PauliSum{N,T}}) where {N,T}
    out = DyadSum(N,T)
    for (dyad, coeff_d) in d
        for (pauli, coeff_p) in p.parent
            new_dyad = dyad*pauli
            sum!(out, new_dyad * coeff_d * coeff_p')
        end   
    end
    return out 
end 
function Base.:*(d::Adjoint{<:Any, DyadSum{N,T}}, p::PauliSum{N,T}) where {N,T}
    out = DyadSum(N,T)
    for (dyad, coeff_d) in d.parent
        for (pauli, coeff_p) in p
            new_dyad = dyad'*pauli
            sum!(out, new_dyad * coeff_d' * coeff_p)
        end   
    end
    return out 
end 

function Base.:*(d::Adjoint{<:Any, DyadSum{N,T}}, p::Adjoint{<:Any, PauliSum{N, T}}) where {N,T}
    out = DyadSum(N,T)
    for (dyad, coeff_d) in d.parent
        for (pauli, coeff_p) in p.parent
            new_dyad = dyad'*pauli
            sum!(out, new_dyad * coeff_d' * coeff_p')
        end   
    end
    return out 
end 



function Base.:*(p::PauliSum{N,T}, d::DyadSum{N,T}) where {N,T}
    out = DyadSum(N,T)
    for (dyad, coeff_d) in d
        for (pauli, coeff_p) in p
            new_dyad = pauli*dyad
            sum!(out, new_dyad * coeff_d * coeff_p)
        end   
    end
    return out 
end 

function Base.:*(p::PauliSum{N,T}, d::Adjoint{<:Any, DyadSum{N,T}}) where {N,T}
    out = DyadSum(N,T)
    for (dyad, coeff_d) in d.parent
        for (pauli, coeff_p) in p
            new_dyad = pauli*dyad'
            sum!(out, new_dyad * coeff_d' * coeff_p)
        end   
    end
    return out 
end 

function Base.:*(p::Adjoint{<:Any, PauliSum{N,T}}, d::DyadSum{N,T}) where {N,T}
    out = DyadSum(N,T)
    for (dyad, coeff_d) in d
        for (pauli, coeff_p) in p.parent
            new_dyad = pauli*dyad
            sum!(out, new_dyad * coeff_d * coeff_p')
        end   
    end
    return out 
end 

function Base.:*(p::Adjoint{<:Any, PauliSum{N,T}}, d::Adjoint{<:Any, DyadSum{N,T}}) where {N,T}
    out = DyadSum(N,T)
    for (dyad, coeff_d) in d.parent
        for (pauli, coeff_p) in p.parent
            new_dyad = pauli*dyad'
            sum!(out, new_dyad * coeff_d' * coeff_p')
        end   
    end
    return out 
end 


promote_to_sum(d::Union{Dyad, DyadBasis}) = DyadSum(d)
promote_to_sum(d::Union{Pauli, PauliBasis}) = PauliSum(d)
promote_to_sum(d::Ket) = KetSum(d)

Singles{N} = Union{Dyad{N}, DyadBasis{N}, Pauli{N}, PauliBasis{N}}
Sums{N,T} = Union{DyadSum{N,T}, PauliSum{N,T}, Adjoint{<:Any, DyadSum{N,T}}, Adjoint{<:Any, PauliSum{N,T}}}

Base.:*(s::Singles{N}, p::Sums{N,T}) where {N,T} = promote_to_sum(s) * p
Base.:*(p::Sums{N,T}, s::Singles{N}) where {N,T} = p * promote_to_sum(s)

