"""
In this representation, the PauliPF string operator is represented as two binary strings, one for x and one for z.
We do not keep track of the phase inside of this type

The format is as follows: 
    
    (i)^θ Z^z1 X^x1 ⊗ Z^z2 X^x2 ⊗ ⋯ ⊗ Z^zN X^xN  
    
Products of operators simply concatonate the left and right strings separately. For example, 

    XYZIY = 11001|01101


To create a Y operator, bits in the same locations in `z` and `x` should be on. 
This means that we have a phase to keep track of because Z^1*X^1 = iY. 
As such, we end up working with operators of the form:
    
    (i)^θ  σ1 ⊗ σ2 ⊗ σ3 ⊗ ⋯ ⊗ σN,

where,

    σ ∈ {X, iY, Z, I}

However, in this case the phase, θ, is fixed to undo the effect of using "iY". 
This means that every `PauliPF` is Hermitian.
The phase angle can be obtained from the function `phase()`
"""
struct PauliPF{N} <: AbstractPauli{N}
    z::Int128
    x::Int128
end


"""
    PauliPF(z::I, x::I) where I<:Integer

TBW
"""
function PauliPF(z::I, x::I, N) where I<:Integer
    # N = maximum(map(i -> ndigits(i, base=2), [x, z]))
    z < 2^N || throw(DimensionMismatch)
    x < 2^N || throw(DimensionMismatch)
    # θ = count_ones(z & x)*3 % 4
    return PauliPF{N}(z, x)
end


"""
    PauliPF(str::String)

Create a `PauliPF` from a string, e.g., 

    a = PauliPF("XXYZIZ")

This is convieniant for manual manipulations, but is not type-stable so will be slow.
"""
function PauliPF(str::String)
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
        # println(i, " ", idx, typeof(idx))
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
    θ = 3*ny%4
    # return PauliPF{N}(z,x), 1im^θ
    return PauliPF{N}(z,x)
end


"""
    PauliPF(N::Integer; X=[], Y=[], Z=[])

constructor for creating PauliPFBoolVec by specifying the qubits where each X, Y, and Z gates exist 
"""
function PauliPF(N::Integer; X=[], Y=[], Z=[])
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
   
    return PauliPF(join(str))
end



"""
    Base.show(io::IO, P::PauliPFMask)

TBW
"""
function Base.show(io::IO, p::PauliPF{N}) where N
    # println(@sprintf "%2i %2iim | %s" real(1im^p.θ) imag(1im^p.θ) string(p)) 
    # println(@sprintf "%s" string(p)) 
    println(@sprintf "%s" string(p)) 
end
"""
    Base.display(p::Pauli)

Display, y = iY
"""
function Base.string(p::PauliPF{N}) where N
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

"""
    random_PauliPF(N)

TBW
"""
function random_PauliPF(N)
    return PauliPF{N}(rand(1:2^N-1),rand(1:2^N-1))
end

"""
    is_hermitian(p::PauliPF)

TBW
"""
function is_hermitian(p::PauliPF)
    return iseven(count_ones(p.x & p.z)) 
end



# @inline get_phase(p::PauliPF) = 1im^phase(p) 
@inline get_phase(p::PauliPF) = 1 
@inline phase(p::PauliPF) = 3*nY(p)%4


"""
    Base.:(==)(p1::Pauli{N}, p2::Pauli{N}) where {N}

Check if they are equal, return true or false
"""
function Base.:(==)(p1::PauliPF, p2::PauliPF) 
    return p1.x == p2.x && p1.z == p2.z
end
function Base.isequal(p1::PauliPF, p2::PauliPF) 
    return p1.x == p2.x && p1.z == p2.z
end


"""
    Base.:(>)(p1::Pauli{N}, p2::Pauli{N}) where {N}

Check if `p1` > `p2`
"""
function Base.:>(p1::PauliPF, p2::PauliPF)
    return p1.z > p2.z || p1.z == p2.z && p1.x > p2.x
end


"""
    Base.:(<)(p1::Pauli{N}, p2::Pauli{N}) where {N}

Check if `p1` < `p2`
"""
function Base.:<(p1::PauliPF, p2::PauliPF)
    return p1.z < p2.z || p1.z == p2.z && p1.x < p2.x
end