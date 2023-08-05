"""
    boson_binary_transformation(nqubits)

TBW
"""
function boson_binary_transformation(nqubits)

    # I = Matrix(Pauli("I"))
    # X = Matrix(Pauli("X"))
    # Y = Matrix(Pauli("Y"))
    # Z = Matrix(Pauli("Z"))

    # σp = (X - 1im*Y) ./ 2
    # σm = (X + 1im*Y) ./ 2
    
    # N = (I + Z) /2
    # M = (I - Z) /2
  
    # out = zeros(2^3, 2^3)
    # out += kron(N,N,σp) .* sqrt(1)
    # out += kron(N,M,σp) .* sqrt(3)
    # out += kron(M,N,σp) .* sqrt(5)
    # out += kron(M,M,σp) .* sqrt(7)
    # out += kron(N,σp,σm) .* sqrt(2)
    # out += kron(M,σp,σm) .* sqrt(6)
    # out += kron(σp,σm,σm) .* sqrt(4)
    # display(out)

    # out = zeros(2^4, 2^4)
    # out += kron(N,N,N,σp) .* sqrt(1)
    # out += kron(N,N,M,σp) .* sqrt(3)
    # out += kron(N,M,N,σp) .* sqrt(5)
    # out += kron(N,M,M,σp) .* sqrt(7)
    # out += kron(M,N,N,σp) .* sqrt(9)
    # out += kron(M,N,M,σp) .* sqrt(11)
    # out += kron(M,M,N,σp) .* sqrt(13)
    # out += kron(M,M,M,σp) .* sqrt(15)
    
    # out += kron(N,N,σp,σm) .* sqrt(2)
    # out += kron(N,M,σp,σm) .* sqrt(6)
    # out += kron(M,N,σp,σm) .* sqrt(10)
    # out += kron(M,M,σp,σm) .* sqrt(14)
    
    # out += kron(N,σp,σm,σm) .* sqrt(4)
    # out += kron(M,σp,σm,σm) .* sqrt(12)
    
    # out += kron(σp,σm,σm,σm) .* sqrt(8)
    # display(out)

    rep1 = []

    for ni in 0:nqubits-1
        stride = 2^(nqubits-ni)

        prodlist = ["I" for i in 1:nqubits]
        prodlist[ni+1] = "σp"
        for nj in ni+2:nqubits
            prodlist[nj] = "σm"
        end
        
        start = 2^(nqubits-ni-1)

        dec_idx = Int(0)

        for i in 1:2^ni
            # str = string(dec_idx, base=2)
            str = bitstring(dec_idx)[64-ni+1:end]
            for (idx,i) in enumerate(str)
                if i == '1'
                    prodlist[idx] = "M" 
                elseif i == '0'
                    prodlist[idx] = "N"
                else
                    throw(ErrorException)
                end 
            end
       
            # println(str)
            push!(rep1, (deepcopy(prodlist), start + dec_idx*stride))
            dec_idx += 1
            
        end
    end

    println(" -------- Intermediate Representation -------- ")
    # println(rep1)
    for i in rep1
        for j in i[1]
            @printf("%2s ",j)
        end
        # println(i[2])
        @printf("√%i\n",i[2])
    end 


    N = ((.5,"I"), (.5, "Z"))
    M = ((.5,"I"), (-.5, "Z"))
    σp = ((.5,"X"), (-.5im, "Y"))
    σm = ((.5,"X"), (.5im, "Y"))
   
    println(" -------- Pauli Representation -------- ")

    rep2 = Dict{Tuple{Int128,Int128},ComplexF64}() 
    # rep2 = Dict{Pauli{nqubits},Float64}() 

    bdag = PauliSum(nqubits)

    for i in rep1
        to_prod = []
        for j in i[1]
            if j == "σp"
                push!(to_prod, σp)
            elseif j == "σm"
                push!(to_prod, σm)
            elseif j == "N"
                push!(to_prod, N)
            elseif j == "M"
                push!(to_prod, M)
            end
        end
        # println()
        # println(i)
        # println(to_prod)
        for prod in Iterators.product(to_prod...)
            coeff = 1
            op = ""
            for i in prod
                coeff *= i[1]
                op *= i[2]
            end
        
            pbs = Pauli(op)
            
            bdag[phasefree(pbs)] = get(bdag, pbs) + get_phase(pbs) * coeff * sqrt(i[2]) 
            
            # if haskey(rep2, key)
            #     rep2[key] += sqrt(i[2]) * coeff * get_phase(pbs) 
            # else
            #     rep2[key] = sqrt(i[2]) * coeff * get_phase(pbs) 
            # end
            # if haskey(rep2, (pbs.z,pbs.x)
            # println(Pauli(op), coeff)
        end
    end
  
    display(bdag)
    return bdag
    for key in sort(rep2)
        @printf("%12.8f %12.8fi %s\n", real(key[2]), imag(key[2]), Pauli(key[1]..., nqubits))
        # println(key[1])
    end
    println(" Number of operators: ", length(rep2))

    rep3 = Dict{Tuple{Int128,Int128},ComplexF64}() 

    println(" Create Product")
    for prod in Iterators.product(keys(rep2), keys(rep2))
        pbs = Pauli(prod[1]..., nqubits) * Pauli(prod[2]..., nqubits)
        # println(pbs)
        coeff::ComplexF64 = rep2[prod[1]]*rep2[prod[2]]'
        key = (pbs.z, pbs.x)
        if haskey(rep3, key)
            rep3[key] += coeff 
            # rep3[key] += coeff * get_phase(pbs) 
        else
            rep3[key] = coeff 
            # rep3[key] = coeff * get_phase(pbs)
        end
    end

    for key in keys(rep3)
        if abs(rep3[key]) < 1e-14
            delete!(rep3,key)
        end
    end
 
    flush(stdout)
    println(" Build b'b matrix")

    mat = zeros(2^nqubits, 2^nqubits)
    for key in sort(rep3)
        pbs = Pauli(key[1]..., nqubits)
        @printf("%12.8f %12.8fi %s\n", real(key[2]), imag(key[2]), pbs)
        mat += key[2]*Matrix(pbs)
        # println(key[1])
    end
    return mat
    
end