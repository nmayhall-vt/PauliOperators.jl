abstract type AbstractPauli{N} end


@inline nqubits(p::AbstractPauli{N}) where N = N
@inline nY(p::AbstractPauli) = count_ones(p.x & p.z)

"""
    Base.Matrix(p::Pauli{N}) where N

Create dense matrix representation 
"""
function Base.Matrix(p::AbstractPauli{N}) where N
    mat = ones(Int8,1,1)
    str = string(p)
    X = [0 1; 1 0]
    y = [0 1; -1 0]
    Y = [0 -1im; 1im 0]
    Z = [1 0; 0 -1]
    I = [1 0; 0 1]
    for i in reverse(1:N)
        # println(i, " ", typeof(str[i]))
        # println(str[i] == "X"[1])
        if str[i] == "X"[1] 
            mat = kron(X,mat)
        elseif str[i] == "y"[1]
            mat = kron(y,mat)
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

    return mat .* get_phase(p)
end


function Base.convert(::Type{Pauli{N}}, p::PauliPF{N}) where {N}
    return Pauli{N}(phase(p), p.z, p.x)
end

function Pauli(p::PauliPF{N}) where {N}
    return Pauli{N}(phase(p), p.z, p.x)
end

"""
    commute(p1::AbstractPauli, p2::AbstractPauli)

Check if they commute, return true or false
"""
function commute(p1::AbstractPauli, p2::AbstractPauli)
    return iseven(count_ones(p1.x & p2.z) - count_ones(p1.z & p2.x)) 
end