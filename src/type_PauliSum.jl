"""
    PauliSum{N, T} = Dict{Tuple{Int128,Int128},T} 

A collection of `Pauli`s, joined by addition.
This uses a `Dict` to store them, however, the specific use cases should probably dictate the container type,
so this will probably be removed.
"""
PauliSum{N, T} = Dict{PauliBasis{N},T} 

PauliSum(N, T) = Dict{PauliBasis{N},T}()

function Base.rand(::Type{PauliSum{N, T}}; n_paulis=2) where {N,T}
    out = PauliSum(N, T)
    for i in 1:n_paulis
        p = rand(Pauli{N})
        out[PauliBasis(p)] = coeff(p) * rand(T)
    end
    return out 
end
function Base.rand(::Type{PauliSum{N}}; n_paulis=2, T=ComplexF64) where {N}
    out = PauliSum(N, T)
    for i in 1:n_paulis
        p = rand(Pauli{N})
        out[PauliBasis(p)] = coeff(p) * rand(T)
    end
    return out 
end

function LinearAlgebra.ishermitian(p::PauliSum{N, T}) where {N,T}
    isherm = true
    for coeff in values(p)
        isherm = isherm && isapprox(imag(coeff), 0, atol=1e-16)
    end
    return isherm
end

function Base.display(ps::PauliSum)
    for (key,val) in ps
        @printf(" %12.8f +%12.8fi %s\n", real(val), imag(val), key)
    end
end

function Base.display(ps::Adjoint{<:Any, PauliSum{N,T}}) where {N,T}
    for (key,val) in ps.parent
        @printf(" %12.8f +%12.8fi %s\n", real(val), -imag(val), key)
    end
end

"""
    Base.Matrix(ps::PauliSum{N}; T=ComplexF64) where N

Create a dense Matrix of type `T` in the standard basis
"""
function Base.Matrix(ps::PauliSum{N, T}) where {N,T}
    out = zeros(T, Int128(2)^N, Int128(2)^N)
    for (op, coeff) in ps
        out .+= Matrix(op) .* coeff 
    end
    return out
end

function LinearAlgebra.tr(p::PauliSum{N, T}) where {N,T}
    return get(p, PauliBasis{N}(0, 0), 0)
end

"""
    Base.:-(ps1::PauliSum, ps2::PauliSum)

Subtract two `PauliSum`s. 
"""
function Base.:-(ps1::PauliSum, ps2::PauliSum)
    out = deepcopy(ps2)
    map!(x->-x, values(out))
    mergewith!(+, out, ps1)
    return out 
end

"""
    Base.:-(ps1::PauliSum, ps2::PauliSum)

Subtract two `PauliSum`s. 
"""
function Base.:-(ps1::PauliSum)
    out = deepcopy(ps1)
    map!(x->-x, values(out))
    return out 
end

Base.adjoint(d::PauliSum{N,T}) where {N,T} = Adjoint(d)
Base.parent(d::Adjoint{<:Any, <:PauliSum}) = d.parent

function Base.Matrix(ps::Adjoint{<:Any, PauliSum{N, T}}) where {N,T}
    out = zeros(T, Int128(2)^N, Int128(2)^N)
    for (op, coeff) in ps.parent
        out .+= Matrix(op) .* adjoint(coeff) 
    end
    return out
end

function Base.size(d::PauliSum{N}) where N
    return (BigInt(2)^N, BigInt(2)^N)
end


function LinearAlgebra.mul!(ps::PauliSum, a::Number)
    map!(x->a*x, values(ps))
    return ps
end


"""
    Base.:*(ps1::PauliSum{N}, ps2::PauliSum{N}) where {N}

Multiply two `PauliSum`s. 
"""
function Base.:*(ps1::PauliSum{N, T}, ps2::PauliSum{N, T}) where {N, T}
    out = PauliSum(N, T)
    for (op1, coeff1) in ps1 
        for (op2, coeff2) in ps2
            prod = Pauli(op1) * Pauli(op2)
            c = coeff(prod)
            prod = PauliBasis(prod)
            if haskey(out, prod)
                out[prod] += c * coeff1 * coeff2
            else
                out[prod] = c * coeff1 * coeff2
            end
        end
    end 
    return out
end
"""
    Base.:*(ps1::Adjoint{<:Any, PauliSum{N, T}}, ps2::PauliSum{N, T}) where {N, T}

Multiply two `PauliSum`s. 
"""
function Base.:*(ps1::Adjoint{<:Any, PauliSum{N, T}}, ps2::PauliSum{N, T}) where {N, T}
    out = PauliSum(N, T)
    for (op1, coeff1) in ps1.parent 
        for (op2, coeff2) in ps2
            prod = Pauli(op1) * Pauli(op2)
            c = coeff(prod)
            prod = PauliBasis(prod)
            if haskey(out, prod)
                out[prod] += c * coeff1' * coeff2
            else
                out[prod] = c * coeff1' * coeff2
            end
        end
    end 
    return out
end
"""
    Base.:*(ps1::PauliSum{N, T}, ps2::Adjoint{<:Any, PauliSum{N, T}}) where {N, T}

Multiply two `PauliSum`s. 
"""
function Base.:*(ps1::PauliSum{N, T}, ps2::Adjoint{<:Any, PauliSum{N, T}}) where {N, T}
    out = PauliSum(N, T)
    for (op1, coeff1) in ps1 
        for (op2, coeff2) in ps2.parent
            prod = Pauli(op1) * Pauli(op2)
            c = coeff(prod)
            prod = PauliBasis(prod)
            if haskey(out, prod)
                out[prod] += c * coeff1 * coeff2'
            else
                out[prod] = c * coeff1 * coeff2'
            end
        end
    end 
    return out
end

"""
    Base.:*(ps1::Adjoint{<:Any, PauliSum{N, T}}, ps2::Adjoint{<:Any, PauliSum{N, T}}) where {N, T}

Multiply two `PauliSum`s. 
"""
function Base.:*(ps1::Adjoint{<:Any, PauliSum{N, T}}, ps2::Adjoint{<:Any, PauliSum{N, T}}) where {N, T}
    out = PauliSum(N, T)
    for (op1, coeff1) in ps1.parent 
        for (op2, coeff2) in ps2.parent
            prod = Pauli(op1) * Pauli(op2)
            c = coeff(prod)
            prod = PauliBasis(prod)
            if haskey(out, prod)
                out[prod] += c * coeff1' * coeff2'
            else
                out[prod] = c * coeff1' * coeff2'
            end
        end
    end 
    return out
end

function Base.:*(ps1::PauliSum{N, T}, a::Number) where {N, T}
    out = deepcopy(ps1)
    mul!(out, a)
    return out
end
Base.:*(a::Number, ps1::PauliSum{N, T}) where {N, T} = ps1 * a

function Base.:*(ps1::Adjoint{<:Any, PauliSum{N, T}}, a::Number) where {N, T}
    out = deepcopy(ps1.parent)
    map!(x->adjoint(x), values(out))
    mul!(out, a)
    return out
end
Base.:*(a::Number, ps1::Adjoint{<:Any, PauliSum{N, T}}) where {N, T} = ps1 * a

Base.getindex(ps::PauliSum, s::String) = ps[PauliBasis(s)]
function Base.getindex(ps::Adjoint{<:Any, PauliSum{N,T}}, a::PauliBasis{N}) where {N,T} 
    return ps.parent[a]'
end

Base.keys(ps::Adjoint{<:Any, PauliSum{N,T}}) where {N,T} = keys(ps.parent)
"""
    Base.sum!(p1::PauliSum{N}, p2::PauliSum{N}) where {N}

Add two `PauliSum`s. 
"""
function Base.sum!(ps1::PauliSum{N}, ps2::PauliSum{N}) where {N}
    mergewith!(+, ps1, ps2)
end
function Base.sum!(ps1::PauliSum{N,T}, ps2::Adjoint{<:Any, PauliSum{N,T}}) where {N,T}
    for (Pauli, coeff) in ps2.parent
        if haskey(ps1, Pauli')
            ps1[Pauli'] += coeff'
        else
            ps1[Pauli'] = coeff'
        end
    end
    return ps1
end
Base.:+(ps2::Adjoint{<:Any, PauliSum{N,T}}, ps1::PauliSum{N,T}) where {N,T} = ps1 + ps2 

"""
    Base.:+(ps1::PauliSum{N}, ps2::PauliSum{N}) where {N}

TBW
"""
function Base.:+(ps1::PauliSum{N}, ps2::PauliSum{N}) where {N}
    out = deepcopy(ps1)
    sum!(out, ps2)
    return out
end
function Base.:+(ps1::PauliSum{N,T}, ps2::Adjoint{<:Any, PauliSum{N,T}}) where {N,T}
    out = deepcopy(ps1)
    sum!(out, ps2)
    return out
end

"""
    otimes(p1::PauliSum{N}, p2::PauliSum{M}) where {N,M}

TBW
"""
function otimes(p1::PauliSum{N,T}, p2::PauliSum{M,T}) where {N,M,T}
    out = PauliSum(N+M, T)
    for (op1,coeff1) in p1
        for (op2,coeff2) in p2
            out[op1 ⊗ op2] = coeff1 * coeff2 
        end
    end
    return out 
end