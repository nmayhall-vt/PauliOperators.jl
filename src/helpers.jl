function get_on_bits(x::T) where T<:Integer
    N = count_ones(x)
    inds = Vector{Int}(undef, N)
    if N == 0
        return inds
    end

    count = 1
    for i in 1:length(bitstring(x))
        if x >> (i-1) & 1 == 1
            inds[count] = i
            count += 1
        end
        count <= N || break
    end
    return inds
end

