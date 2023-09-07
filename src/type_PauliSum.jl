using LinearAlgebra
using BlockDavidson 

"""
   ::Dict{Pauli{N},ComplexF64}

A collection of `Pauli`s, joined by addition.
This uses a `Dict` to store them, however, the specific use cases should probably dictate the container type,
so this will probably be removed.
"""
PauliSum{N} = Dict{Pauli{N},ComplexF64}


"""
    PauliSum(N)

TBW
"""
function PauliSum(N)
    return PauliSum{N}(Dict{Pauli{N},ComplexF64}())
end


"""
    Base.display(ps::PauliSum)

TBW
"""
function Base.display(ps::PauliSum)
    for (key,val) in ps
        @printf(" %12.8f +%12.8fi θ:%1i %s\n", real(val), imag(val), key.θ, key)
    end
end

Base.get(ps::PauliSum{N}, p::Pauli{N}) where N = get(ps, phasefree(p), zero(ComplexF64))
Base.getindex(ps::PauliSum{N}, p::Pauli{N}) where N = ps[phasefree(p)]
Base.setindex!(ps::PauliSum{N}, v, p::Pauli{N}) where N = ps[phasefree(p)] = v*get_phase(p)

"""
    Base.+(p1::PauliSum{N}, p2::PauliSum{N}) where {N}

Add two `PauliSum`s. 
"""
function Base.:+(ps1::PauliSum{N}, ps2::PauliSum{N}) where {N}
    return PauliSum{N}(mergewith(+, ps1, ps2))
end

"""
    Base.sum!(p1::PauliSum{N}, p2::PauliSum{N}) where {N}

Add two `PauliSum`s. 
"""
function Base.sum!(ps1::PauliSum{N}, ps2::PauliSum{N}) where {N}
    mergewith!(+, ps1, ps2)
end

"""
    Base.-(p1::PauliSum{N}, p2::PauliSum{N}) where {N}

Subtract two `PauliSum`s. 
"""
function Base.:-(ps1::PauliSum{N}, ps2::PauliSum{N}) where {N}
    return PauliSum{N}(mergewith(-, ps1, ps2))
end

Base.length(ps::PauliSum) = length(ps)
"""
    LinearAlgebra.adjoint!(ps::PauliSum{N}) where N

TBW
"""
function LinearAlgebra.adjoint!(ps::PauliSum{N}) where N 
    for (key, val) in ps
        # ps[key] = adjoint(val)  
        ps[key] = adjoint(val) * (-1)^count_ones(key.z & key.x) 
    end
end
"""
    LinearAlgebra.adjoint(ps::PauliSum{N}) where N

TBW
"""
function LinearAlgebra.adjoint(ps::PauliSum{N}) where N 
    out = deepcopy(ps)
    adjoint!(out)    
    return out
end

#
# lazy adjoint views aren't yet working. 
#
# Base.adjoint(ps::PauliSum{N}) where N = Adjoint{ComplexF64,PauliSum{N}}(ps)
# Base.size(aps::Adjoint{ComplexF64, PauliSum{N}}) where N = (length(aps),)
# function Base.adjoint(ps::PauliSum)
#     return LinearAlgebra.Adjoint
# end

"""
    Base.-(p1::PauliSum{N}, p2::PauliSum{N}) where {N}

Multiply two `PauliSum`s. 
"""
function Base.:*(ps1::PauliSum{N}, ps2::PauliSum{N}) where {N}
    out = PauliSum(N)
    for (op1, coeff1) in ps1 
        for (op2, coeff2) in ps2
            prod = op1 * op2
            out[phasefree(prod)] = get(out, prod) + get_phase(prod)*coeff1*coeff2 
        end
    end 
    return out
end

"""
    Base.:*(ps::PauliSum{N}, a::Number) where {N}

TBW
"""
function Base.:*(ps::PauliSum{N}, a::Number) where {N}
    out = deepcopy(ps) 
    mul!(out,a)
    return out
end
Base.:*(a::Number, ps::PauliSum{N}) where {N} = ps*a

"""
    LinearAlgebra.mul!(ps::PauliSum, a::Number)

TBW
"""
function LinearAlgebra.mul!(ps::PauliSum, a::Number)
    for (op, coeff) in ps 
        ps[op] = coeff * a
    end 
end

"""
    Base.Matrix(ps::PauliSum{N}; T=ComplexF64) where N

Create a dense Matrix of type `T`
"""
function Base.Matrix(ps::PauliSum{N}; T=ComplexF64) where N
    out = zeros(T, 2^N, 2^N)
    for (op, coeff) in ps
        out .+= Matrix(op) .* coeff
    end
    return out
end

"""
    clip!(ps::PauliSum; thresh=1e-16)

Delete Pauli's with coeffs smaller than thresh
"""
function clip!(ps::PauliSum; thresh=1e-16)
    to_delete = []
    for (op,coeff) in ps
        if abs(coeff) < thresh
            push!(to_delete, op)
        end
    end
    for k in to_delete
        delete!(ps, k)
    end
end

"""
    Base.:≈(p1::PauliSum{N}, p2::PauliSum{N}) where {N}

TBW
"""
function Base.:≈(p1::PauliSum{N}, p2::PauliSum{N}) where {N}
    for (op, coeff) in p1
        get(p2, op) .≈ coeff || return false
    end
    for (op, coeff) in p2
        get(p1, op) .≈ coeff || return false
    end
    return true
end


function BlockDavidson.LinOpMat(ps::PauliSum{N}; issymmetric=false) where N
    function matvec(v)
        return ps * v
    end
    return LinOpMat{ComplexF64}(matvec, 2^N, issymmetric) 
end


# function LinearMaps.LinearMap(ps::PauliSum{N}) where N

#     function matvec(v)
#         return ps * v
#     end
#     return LinearMap 
# end

function nqubits(p::PauliSum{N}) where N
    return N
end