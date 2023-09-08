## 
# These things are not yet defined, but might be considered in the future

"""
    θ::UInt8
    z::NTuple{R,Int32}
    x::NTuple{R,Int32}
"""
struct PauliNew{N,R}
    θ::UInt8
    z::NTuple{R,Int32}
    x::NTuple{R,Int32}
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
    return PauliNew{N,M}(θ, ntuple(zints->zints, M), ntuple(xints->xints, M))
end