
function Base.Matrix(d::Dyad{N}) where N
    mat = zeros(Bool, size(d))
    mat[d.ket.v+1, d.bra.v+1] = true
    return mat 
end

function Base.Matrix(d::ScaledDyad{N,T}) where {N,T}
    mat = zeros(T, size(d))
    mat[d.dyad.ket.v+1, d.dyad.bra.v+1] = d.coeff 
    return mat 
end

function Base.Matrix(d::DyadSum{N,T}) where {N,T}
    mat = zeros(T, size(d))
    for (dyad, coeff) in d
        mat[dyad.ket.v+1, dyad.bra.v+1] = coeff 
    end
    return mat 
end

function Base.Matrix(p::FixedPhasePauli{N}) where N
    mat = ones(Int8,1,1)
    str = string(p)
    X = [0 1; 1 0]
    y = [0 1; -1 0]
    Y = [0 -1im; 1im 0]
    Z = [1 0; 0 -1]
    I = [1 0; 0 1]
    # for i in reverse(1:N)
    for i in 1:N
        if str[i] == "X"[1] 
            mat = kron(X,mat)
        elseif str[i] == "y"[1]
            mat = kron(y,mat)
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

function Base.Matrix(d::AbstractState{N}) where N
    mat = zeros(Bool, size(d))
    mat[d.v+1] = true
    return mat 
end

"""
    Base.Matrix(p::Pauli)

Create dense matrix representation in standard basis 
"""
Base.Matrix(p::Pauli) = Matrix(p.pauli) .* get_phase(p)

"""
    Base.Matrix(ps::PauliSum{N}; T=ComplexF64) where N

Create a dense Matrix of type `T` in the standard basis
"""
function Base.Matrix(ps::PauliSum{N}; T=ComplexF64) where N
    out = zeros(T, Int128(2)^N, Int128(2)^N)
    for (op, coeff) in ps.ops
        out .+= Matrix(op) .* coeff 
    end
    return out
end

"""
    Base.Matrix(spv::Vector{ScaledPauli{N}}) where {T,N}

Create dense matrix representation in standard basis 
"""
function Base.Matrix(spv::Vector{ScaledPauli{N}}; T=ComplexF64) where {N}
    out = zeros(T, Int128(2)^N, Int128(2)^N)
    for spvi in spv 
        out .+= Matrix(spvi.pauli) .* spvi.coeff
    end
    return out
end


"""
    Base.Matrix(p::Pauli{N}) where N

Create dense matrix representation in standard basis 
"""
function Base.Matrix(p::ScaledPauli{N}) where {N}
    return Matrix(p.pauli) .* p.coeff 
end
