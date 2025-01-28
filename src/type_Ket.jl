"""
An occupation number vector, up to 128 qubits
"""
struct Ket{N} 
    v::Int128
end

struct Bra{N}
    v::Int128
end

KetSum{N, T} = Dict{Ket{N}, T}

"""
    Ket(vec::Vector{T}) where T<:Union{Bool, Integer}

TBW
"""
function Ket(vec::Vector{T}) where T<:Union{Bool, Integer}
    two = Int128(2)
    v = Int128(0)

    for i in 1:length(vec)
        if vec[i] == 1
            v |= two^(i-1)
        end
    end
    return Ket{length(vec)}(v)
end
function Bra(vec::Vector{T}) where T<:Union{Bool, Integer}
    two = Int128(2)
    v = Int128(0)

    for i in 1:length(vec)
        if vec[i] == 1
            v |= two^(i-1)
        end
    end
    return Bra{length(vec)}(v)
end

"""
    Ket(N::Integer, v::Integer)

TBW
"""
function Ket(N::Integer, v::Integer)
    for i in N+1:128
        v &= ~(Int128(2)^(i-1))    
    end
    return Ket{N}(Int128(v))
end
function Bra(N::Integer, v::Integer)
    for i in N+1:128
        v &= ~(Int128(2)^(i-1))    
    end
    return Bra{N}(Int128(v))
end


function KetSum(N; T=Float64)
    return Dict{Ket{N}, T}()
end
    
function Base.size(d::Ket{N}) where N
    return (BigInt(2)^N, 1)
end

function Base.size(d::Bra{N}) where N
    return (1, BigInt(2)^N)
end

Base.adjoint(d::Ket{N}) where N = Bra{N}(d.v)
Base.adjoint(d::Bra{N}) where N = Ket{N}(d.v)


Base.rand(T::Type{Ket{N}}) where N = T(rand(0:Int128(2)^N-1))
Base.rand(T::Type{Bra{N}}) where N = T(rand(0:Int128(2)^N-1))




"""
    Base.show(io::IO, P::Pauli{N}) where N

TBW
"""
function Base.show(io::IO, P::Union{Ket,Bra})
    print(io, string(P))
end

function Base.string(p::Ket{N}) where N
    out = [0 for i in 1:128]
    for i in get_on_bits(p.v)
        out[i] = 1
    end
    return "|"*join(out[1:N])*">"
end
function Base.string(p::Bra{N}) where N
    out = [0 for i in 1:128]
    for i in get_on_bits(p.v)
        out[i] = 1
    end
    return "<"*join(out[1:N])*"|"
end


"""
    Base.:+(p::Ket{N}, q::Ket{N}) where N

Add two `Ket`'s together, return a `KetSum`
"""
function Base.:+(p::Ket{N}, q::Ket{N}) where N
    if p == q
        return KetSum{N, ComplexF64}(p => 2) 
    else 
        return KetSum{N, ComplexF64}(p => 1, q => 1)
    end
end
"""
    Base.:+(p::Bra{N}, q::Bra{N}) where N

Add two `Ket`'s together, return a `KetSum`
"""
function Base.:+(p::Bra{N}, q::Bra{N}) where N
    if p == q
        return KetSum{N, ComplexF64}(p' => 2)'
    else 
        return KetSum{N, ComplexF64}(p' => 1, q' => 1)'
    end
end



function index(k::Union{Ket{N}, Bra{N}}) where N
    return k.v+1
end

"""
    Base.Vector(k::Union{Ket{N}, Bra{N}}; T=Int64) where N

Create dense vector representation in standard basis 
"""
function Base.Vector(k::Union{Ket{N}, Bra{N}}; T=Int64) where N
    vec = zeros(T,Int128(2)^N)
    vec[index(k)] = T(1) 
    return vec 
end


function Base.iterate(::Type{Ket{N}}, state = 1) where N
    state > 4^N && return
    return Ket{N}(state-1), state+1 
end
