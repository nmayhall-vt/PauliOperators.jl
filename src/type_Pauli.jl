"""
In this representation, the Pauli string operator is represented as two binary strings, one for x and one for z.

The format is as follows: 
    
    (i)^\theta Z^z1 X^x1 ⊗ Z^z2 X^x2 ⊗ ⋯ ⊗ Z^zN X^xN  
    
Products of operators simply concatonate the left and right strings separately. For example, 

    XYZIY = 11001|01101


To create a Y operator, bits in the same locations in `z` and `x` should be on. 
This means that we have a phase to keep track of because Z^1*X^1 = iY. 
As such, we end up working with operators of the form:
    
    (i)^\theta σ1 ⊗ σ2 ⊗ σ3 ⊗ ⋯ ⊗ σN,

where,

    σ ∈ {X, iY, Z, I}

"""
struct Pauli{N} <: AbstractPauli{N}
    θ::UInt8
    pauli::FixedPhasePauli{N}
end


"""
    Pauli(z::I, x::I) where I<:Integer

TBW
"""
function Pauli(z::I, x::I, N) where I<:Integer
    # N = maximum(map(i -> ndigits(i, base=2), [x, z]))
    z < Int128(2)^N || throw(DimensionMismatch)
    x < Int128(2)^N || throw(DimensionMismatch)
    θ = count_ones(z & x)*3 % 4
    return Pauli{N}(θ, FixedPhasePauli{N}(z, x))
end


"""
    Pauli(str::String)

Create a `Pauli` from a string, e.g., 

    a = Pauli("XXYZIZ")

This is convieniant for manual manipulations, but is not type-stable so will be slow.
"""
function Pauli(str::String)
    p = FixedPhasePauli(str)
    return Pauli{nqubits(p)}(phase(p), p) 
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
    PauliNew(N::Integer; X=[], Y=[], Z=[])

constructor for creating PauliBoolVec by specifying the qubits where each X, Y, and Z gates exist 
"""
function PauliNew(N::Integer; X=[], Y=[], Z=[])
    M = (N-1)÷32+1
    one = Int32(1)
    two = Int32(2)


    # println(M)
    for i in X
        i ∉ Y || throw(DimensionMismatch)
        i ∉ Z || throw(DimensionMismatch)
    end
    for i in Y
        i ∉ Z || throw(DimensionMismatch)
    end
   
    zints = [Int32(0) for i in 1:M]
    xints = [Int32(0) for i in 1:M]
    for i in X
        register = (i-1)÷32+1
        idx = i%32+1
        xints[register] |= two^(idx-one)
        # println(register, " ", index)
    end
    for i in Y
        register = (i-1)÷32+1
        idx = i%32+1
        zints[register] |= two^(idx-one)
        xints[register] |= two^(idx-one)
    end
    for i in Z
        register = (i-1)÷32+1
        idx = i%32+1
        zints[register] |= two^(idx-one)
    end
    
    θ = 3*length(Y)%4 
    return PauliNew{N,M}(θ, FixedPhasePauli{N}(ntuple(zints->zints, M), ntuple(xints->xints, M)))
end


"""
    Base.show(io::IO, P::PauliMask)

TBW
"""
function Base.show(io::IO, p::Pauli{N}) where N
    println(@sprintf "%2i %2iim | %s" real(1im^p.θ) imag(1im^p.θ) string(p)) 
end

"""
    Base.display(p::Pauli)

Display, y = iY
"""
Base.string(p::Pauli{N}) where N = string(p.pauli)

"""
    Base.Matrix(p::Pauli)

TBW
"""
Base.Matrix(p::Pauli) = Matrix(p.pauli) .* get_phase(p)

"""
    random_Pauli(N)

TBW
"""
function random_Pauli(N)
    return Pauli{N}(rand(0:3), random_FixedPhasePauli(N))
end

"""
    is_hermitian(p::Pauli)

TBW
"""
function is_hermitian(p::Pauli)
    return ~(iseven(p.θ) ⊻ is_hermitian(p.pauli))
end



"""
    Base.:-(p::Pauli{N}) where {N}

TBW
"""
function Base.:-(p::Pauli{N}) where {N}
    return rotate_phase(p,2) 
end


Base.adjoint(p::Pauli) = is_hermitian(p) ? p : -p

LinearAlgebra.tr(p::Pauli) = get_phase(p) * tr(p.pauli)