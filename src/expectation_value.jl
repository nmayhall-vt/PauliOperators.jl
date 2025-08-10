
function expectation_value(p::PauliBasis{N}, ket::Ket{N}) where N
    is_diagonal(p) || return 0.0

    count_ones(p.z & ket.v) % 2 == 0 || return -1
    return 1
end



function expectation_value(p::Union{PauliBasis{N}, Pauli{N}}, d::Union{Dyad{N}, DyadBasis{N}}) where N
    sgn = count_ones(p.z & d.bra.v)  # sgn <j| = <j| z 
    val = d.ket.v ⊻ d.bra.v == p.x # <j|x|i>
    return (-1)^sgn * val * coeff(p) * coeff(d) * 1im^symplectic_phase(p)
end

function expectation_value(p::PauliSum{N,T}, d::Union{Dyad{N}, DyadBasis{N}}) where {N,T}
    eval = zero(T)
    for (pi,ci) in p
        eval += expectation_value(pi, d) * ci
    end
    return eval 
end

function expectation_value(p::Union{PauliBasis{N}, Pauli{N}}, d::DyadSum{N,T}) where {N,T}
    eval = zero(T)
    for (di,ci) in d
        eval += expectation_value(p, di) * ci
    end
    return eval 
end

function expectation_value(p::PauliSum{N,T}, d::DyadSum{N,T}) where {N,T}
    eval = zero(T)
    for (pi,ci) in p
        for (dj,cj) in d
            eval += expectation_value(pi, dj) * ci * cj
        end
    end
    return eval 
end
