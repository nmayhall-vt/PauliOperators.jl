
"""
    Base.:*(p1::PauliBitString{N}, p2::PauliBitString{N}) where {N}

Multiply two `PauliBitString`'s together
"""
function Base.:*(p1::PauliBitString{N}, p2::PauliBitString{N}) where {N}
    x = p1.x ⊻ p2.x
    z = p1.z ⊻ p2.z
    θ = (p1.θ + p2.θ ) % 4
    θ = (θ + 2*count_ones(p1.x & p2.z)) % 4
    return PauliBitString{N}(θ,z,x)
end


"""
    Base.:(==)(p1::PauliBitString{N}, p2::PauliBitString{N}) where {N}

Check if they are equal, return true or false
"""
function Base.:(==)(p1::PauliBitString{N}, p2::PauliBitString{N}) where {N}
    return p1.x == p2.x && p1.z == p2.z && p1.θ == p2.θ
end



"""
    get_phase(p::PauliBitString)

Return the phase of the `PauliBitString`, i^θ
"""
function get_phase(p::PauliBitString)
    return 1im^p.θ
end

"""
    rotate_phase(p::PauliBitString{N}, θ::Integer) where N

Rotate phase in units of π/2. In otherwords, multiply the phase by i^θ.
E.g., mutliplication by -1 is obtained with θ=2.
"""
function rotate_phase(p::PauliBitString{N}, θ::Integer) where N
    return PauliBitString{N}((p.θ + θ)%4, p.z, p.x)
end


"""
    negate(p::PauliBitString)

Multiply `p` by -1
"""
function negate(p::PauliBitString{N}) where N
    return rotate_phase(p,2) 
end

"""
    is_diagonal(p::PauliBitString)

Check if operator is diagonal in the computational (z) basis. E.g., does this operator consist of only I and/or Z?
"""
function is_diagonal(p::PauliBitString)
    return count_ones(p.x) == 0
end

"""
    Base.show(io::IO, P::PauliMask)

TBW
"""
function Base.show(io::IO, p::PauliBitString{N}) where N
    # print(io, @sprintf "Pstring(P))
    # println(io, 1im^p.θ,"|", string(p)) 
    println(@sprintf "%2i %2iim | %s" real(1im^p.θ) imag(1im^p.θ) string(p)) 
end


"""
    Base.display(p::PauliBitString)

Display, y = iY
"""
function Base.string(p::PauliBitString{N}) where N
    Iloc = get_on_bits(p.x ⊽ p.z)
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
    random_PauliBitString(N)

TBW
"""
function random_PauliBitString(N)
    return PauliBitString{N}(rand(0:3), rand(Int128),rand(Int128))
end


"""
    commute(p1::PauliBitString{N}, p2::PauliBitString{N}) where {N}

Check if they commute, return true or false
"""
function commute(p1::PauliBitString{N}, p2::PauliBitString{N}) where {N}
    return iseven(count_ones(p1.x & p2.z) - count_ones(p1.z & p2.x)) 
end


"""
    expectation_value_sign(p::PauliBitString{N}, ket::Vector{Bool}) where N

compute expectation value of PauliBitString `o` for a product state `ket`
"""
function expectation_value_sign(p::PauliBitString{N}, ket::KetBitString{N}) where N
    is_diagonal(p) || return 0.0
    
    println(p)

    count_ones(p.z & ket.v) % 2 == 0 || return -(1im)^p.θ
    return (1im)^p.θ 
end
