
"""
    z::Int128
    x::Int128

A positive, Hermitian Pauli, used as a basis for more general `Pauli`'s (which can have a complex phase).
These are primarily used to provide a basis for linear combinations of Paulis, e.g., `PauliSum`'s.
    
    PauliBasis{N}(z,x)  =  i^θs ⋅ z₁...|x₁... 
                            =  P₁⊗...⊗Pₙ

Phase definitions:
- `symplectic_phase`: `θs` - phase needed to cancel the phase arising from the ZX factorized form: `θs = θ-θg`
"""
struct PauliBasis{N} 
    z::Int128
    x::Int128
end

LinearAlgebra.ishermitian(p::PauliBasis) = true
coeff(p::PauliBasis) = 1

@inline symplectic_phase(p::PauliBasis) = (4-count_ones(p.z & p.x)%4)%4

function PauliBasis(str::String)
    for i in str
        i in ['I', 'Z', 'X', 'Y'] || error("Bad string: ", str)
    end

    x = Int128(0)
    z = Int128(0)
    ny = 0 
    N = length(str)
    idx = Int128(1)
    two = Int128(2)
    one = Int128(1)

    for i in str
        if i in ['X', 'Y']
            x |= two^(idx-one)
            if i == 'Y'
                ny += 1
            end
        end
        if i in ['Z', 'Y']
            z |= two^(idx-one)
        end
        idx += 1
    end
    return PauliBasis{N}(z, x)
end



"""
    Base.Matrix(p::PauliBasis{N}) where N

Build dense matrix representation in standard basis
"""
function Base.Matrix(p::PauliBasis{N}) where N
    mat = ones(Int8,1,1)
    str = string(p)
    X = [0 1; 1 0]
    Y = [0 -1im; 1im 0]
    Z = [1 0; 0 -1]
    I = [1 0; 0 1]
    # for i in reverse(1:N)
    for i in 1:N
        if str[i] == "X"[1] 
            mat = kron(X,mat)
        elseif str[i] == "Y"[1]
            mat = kron(Y,mat)
        elseif str[i] == "Z"[1]
            mat = kron(Z,mat)
        elseif str[i] == "I"[1]
            mat = kron(I,mat)
        else
            throw(ErrorException)
        end
    end

    return mat
end

PauliBasis(p::PauliBasis) = p

"""
    Base.string(p::Pauli{N}) where N

Display, y = iY
"""
function Base.string(p::PauliBasis{N}) where N
    yloc = get_on_bits(p.x & p.z)
    Xloc = get_on_bits(p.x & ~p.z)
    Zloc = get_on_bits(p.z & ~p.x)
    out = ["I" for i in 1:128]

    for i in Xloc
        out[i] = "X"
    end
    for i in yloc
        out[i] = "Y"
    end
    for i in Zloc
        out[i] = "Z"
    end
    return join(out[1:N])
end

function Base.rand(T::Type{PauliBasis{N}}) where N
    return PauliBasis{N}(rand(0:Int128(2)^N-1), rand(0:Int128(2)^N-1))
end

Base.display(p::PauliBasis) = println(string(p))

function otimes(p1::PauliBasis{N}, p2::PauliBasis{M}) where {N,M} 
    PauliBasis{N+M}(p1.z | p2.z << N, p1.x | p2.x << N)
end

Base.:*(p1::PauliBasis, p2::PauliBasis) = Pauli(p1) * Pauli(p2)
Base.:*(p1::PauliBasis{N}, a::Number) where N = Pauli(p1)*a
Base.:*(a::Number, p1::PauliBasis{N}) where N = Pauli(p1)*a

Base.adjoint(p::PauliBasis) = p

function Base.iterate(::Type{PauliBasis{N}}, state = 1) where N
    state > 4^N && return
    next = CartesianIndices((2^N,2^N))[state]
    return PauliBasis{N}(next[1]-1, next[2]-1), state+1 
end
 
@inline commute(p1::PauliBasis, p2::PauliBasis) = iseven(count_ones(p1.x & p2.z) - count_ones(p1.z & p2.x)) 