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
struct Pauli{N} <: Integer
    θ::UInt8
    z::Int128
    x::Int128
end

"""
    Pauli(z::I, x::I) where I<:Integer

TBW
"""
function Pauli(z::I, x::I, N) where I<:Integer
    # N = maximum(map(i -> ndigits(i, base=2), [x, z]))
    z < 2^N || throw(DimensionMismatch)
    x < 2^N || throw(DimensionMismatch)
    θ = count_ones(z & x)*3 % 4
    return Pauli{N}(θ, z, x)
end

"""
    Pauli(str::String)

TBW
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
    return Pauli{N}(θ, z,x) 
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

function phasefree(p::Pauli)
    return rotate_phase(p, -p.θ)
end