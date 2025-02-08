"""
    Pauli{N} 

is our basic type for representing Pauli operators acting on `N`.
Assume we want to represent a Pauli string of the following form:

    σ1 ⊗ σ2 ⊗ σ3 ⊗ ⋯ ⊗ σN,

where, `σ ∈ {X, Y, Z, I}`.
To do this efficiently, we use the symplectic representation of the Pauli group, 
where we factor each Pauli into a product of X and Z operators: 
    
    σ = i^(3*(z+x)%2) Zᶻ Xˣ,
    
with z,x ∈ {0,1}. The phase factor comes from the fact that `Z*X = iY`.
In this representation, any tensor product of Pauli's is represented as two binary strings, one for x and one for z, along with the associated phase accumulated from each site.
The format is as follows: 
    
    i^θ   Z^z₁ ⋅ X^x₁ ⊗ Z^z₂ ⋅ X^x₂ ⊗ ⋯ ⊗ Z^zₙ ⋅ X^xₙ  
    
Products of operators simply concatonate the left and right strings separately. For example, 
To create a Y operator, bits in the same locations in `z` and `x` should be on. 
    
    XYZIy = 11001|01101     where y = iY

Since we get a factor of `i` each time we create a Y operator, we need to keep track of this to cancel the  phase `θs`, arising from the ZX factorization.
    
    P₁⊗...⊗Pₙ = i^θs ⋅ z₁...|x₁...  where Pᵢ ∈ {I,X,Y,Z}.

similarly, 

    z₁...|x₁... = i^-θs ⋅ P₁⊗...⊗Pₙ 

We use `θs` to denote the phase needed to make the Pauli operator Hermitian and positive, and we refer to this as the `symplectic_phase`, since it arises solely from the symplectic representation of the Pauli.
However, this is not the only phase we need to worry about. Since various phases accumulate during Pauli multiplication, we allow a given `Pauli` to have an arbitrary global phase, `θg`, so that the `Pauli` type can be closed under multiplication. As such, our `Pauli` phases are defined according to the following:

    Pauli{N}(s,z,x)  =  s ⋅ z₁...|x₁... 
                     =  s ⋅ i^-θs ⋅ P₁⊗...⊗Pₙ
                     =  coeff ⋅ P₁⊗...⊗Pₙ

    PauliBasis{N}(z,x)  =  i^θs ⋅ z₁...|x₁... 
                            =  P₁⊗...⊗Pₙ

Phase definitions:
- `symplectic_phase`: `θs` - phase needed to cancel the phase arising from the ZX factorized form: `θs = θ-θg`

Since we need to keep track of a phase for a Pauli, we might as well let it become a general scalar value for broader use. As such, `Pauli.s` is a arbitrary complex number.
"""
struct Pauli{N} 
    s::ComplexF64
    z::Int128
    x::Int128
end

PauliTypes{N} = Union{Pauli{N}, PauliBasis{N}}

"""
    coeff(p::Pauli)

Return the coefficient from the product of the scalar times the inverse symplectic_phase
"""
# coeff(p::Pauli) = 1im^Float64((4 - symplectic_phase(p) )%4)
# symplectic_phase(p::Pauli) = (4-count_ones(p.z & p.x)%4)%4
@inline coeff(p::Pauli) = p.s * 1im^inv_symplectic_phase(p)
@inline inv_symplectic_phase(p::Pauli) = (4-symplectic_phase(p)%4)
@inline symplectic_phase(p::Pauli) = (4-count_ones(p.z & p.x)%4)%4

function Pauli(p::PauliBasis{N}) where N
    return Pauli{N}(1im^symplectic_phase(p), p.z, p.x)
end

"""
    Pauli(z::I, x::I) where I<:Integer

TBW
"""
function Pauli(z::I, x::I, N) where I<:Integer
    # N = maximum(map(i -> ndigits(i, base=2), [x, z]))
    z < Int128(2)^N || throw(DimensionMismatch)
    x < Int128(2)^N || throw(DimensionMismatch)
    # θ = get_hermitian_phase(z,x)
    return Pauli{N}(1, z, x)
end


"""
    Pauli(str::String)

Create a `Pauli` from a string, e.g., 

    a = Pauli("XXYZIZ")

This is convieniant for manual manipulations, but is not type-stable so will be slow.
"""
function Pauli(str::String)
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
    θ = 4-ny%4
    return Pauli{N}(1im^θ, z, x)
end




"""
    Pauli(N::Integer; X=[], Y=[], Z=[])

constructor for creating PauliBoolVec by specifying the qubits where each X, Y, and Z gates exist 
"""
function Pauli(N::Integer; X=[], Y=[], Z=[])
    for i in X
        i ∉ Y || throw(DimensionMismatch)
        i ∉ Z || throw(DimensionMismatch)
    end
    for i in Y
        i ∉ Z || throw(DimensionMismatch)
    end
    
    str = ["I" for i in 1:N]
    for i in X
        str[i] = "X"
    end
    for i in Y
        str[i] = "Y"
    end
    for i in Z
        str[i] = "Z"
    end
   
    # print(str[1:N])
    return Pauli(join(str))
end   




"""
    Base.string(p::Pauli{N}) where N

Display, y = iY
"""
function Base.string(p::Pauli{N}) where N
    yloc = get_on_bits(p.x & p.z)
    Xloc = get_on_bits(p.x & ~p.z)
    Zloc = get_on_bits(p.z & ~p.x)
    out = ["I" for i in 1:128]

    for i in Xloc
        out[i] = "X"
    end
    for i in yloc
        out[i] = "y"
    end
    for i in Zloc
        out[i] = "Z"
    end
    return join(out[1:N])
end

"""
    Base.display(p::Pauli{N}) where N

TBW
"""
function Base.display(p::Pauli{N}) where N
    @printf "% .4f % .4fim | %s\n" real(p.s) imag(p.s) string(p) 
end


"""
    rand(Pauli{N})

TBW
"""
function Base.rand(T::Type{Pauli{N}}) where N
    return Pauli{N}(randn(ComplexF64), rand(0:Int128(2)^N-1), rand(0:Int128(2)^N-1))
end

function nY(p::Pauli)
    return count_ones(p.x & p.z)
end

"""
    ishermitian(p::Pauli)

TBW
"""
function LinearAlgebra.ishermitian(p::Pauli; thresh=1e-16)
    return abs(imag(coeff(p)))<thresh
end


"""
    Base.Matrix(p::Pauli{N}) where N

Build dense matrix representation in standard basis
"""
Base.Matrix(p::Pauli) = Matrix(PauliBasis(p)) * coeff(p)


"""
    Base.:-(p::Pauli{N}) where {N}

TBW
"""
function Base.:-(p::Pauli{N}) where {N}
    return Pauli{N}(-p.s, p.z, p.x) 
end


"""
    Base.adjoint(p::Pauli)


    Pauli{N}(s,z,x)  =  s ⋅ z₁...|x₁... 
                     =  s ⋅ i^-θs ⋅ P₁⊗...⊗Pₙ
                     =  coeff ⋅ P₁⊗...⊗Pₙ

Since the PauliBasis is Hermitian, we have that
    Pauli' = coeff' ⋅ P₁⊗...⊗Pₙ
"""
Base.adjoint(p::Pauli{N}) where N = Pauli{N}(coeff(p)'*1im^symplectic_phase(p), p.z, p.x)
# Base.adjoint(p::Pauli{N}) where N = Pauli{N}(coeff(p)'*1im^symplectic_phase(p), p.z, p.x)

function LinearAlgebra.tr(p::Union{Pauli{N}, PauliBasis{N}}) where N
    return coeff(p) * ((p.z == 0) && (p.x == 0)) * 2^N
end

"""
    Base.:*(p1::Pauli{N}, p2::Pauli{N}) where {N}

Multiply two `Pauli`'s together
"""
function Base.:*(p1::Pauli{N}, p2::Pauli{N}) where {N}
    x = p1.x ⊻ p2.x
    z = p1.z ⊻ p2.z
    s = p1.s * p2.s * 1im^(2*count_ones(p1.x & p2.z) % 4)
    return Pauli{N}(s, z, x)
end

Base.:*(p::Pauli{N}, s::Number) where N = Pauli{N}(p.s * s, p.z, p.x)
Base.:*(s::Number, p::Pauli{N}) where N = p*s 

"""
    Base.:+(p::Pauli{N}, q::Pauli{N}) where N

Add two `Pauli`'s together, return a `PauliSum`
"""
function Base.:+(p::PauliTypes{N}, q::PauliTypes{N}) where N
    if PauliBasis(p) == PauliBasis(q)
        return PauliSum{N, ComplexF64}(PauliBasis(p) => coeff(p)+coeff(q))
    else 
        return PauliSum{N, ComplexF64}(PauliBasis(p) => coeff(p), PauliBasis(q) => coeff(q))
    end
end

"""
    otimes(p1::Pauli{N}, p2::Pauli{M}) where {N,M}

TBW
"""
function otimes(p1::Pauli{N}, p2::Pauli{M}) where {N,M} 
    Pauli{N+M}(p1.s * p2.s, p1.z | p2.z << N, p1.x | p2.x << N)
end

"""
    osum(p1::Pauli{N}, p2::Pauli{M}) where {N,M}

Returns the direct sum of two Paulis
"""
function osum(p1::Pauli{N}, p2::Pauli{M}) where {N,M}
    return p1 ⊗ Pauli(M) + Pauli(N) ⊗ p2 
end

function Base.iterate(::Type{Pauli{N}}, state = 1) where N
    state > 4^N && return
    next = CartesianIndices((2^N,2^N))[state]
    bp = PauliBasis{N}(next[1]-1, next[2]-1)
    return Pauli(bp), state+1 
end