function boson_binary_transformation(nqubits)

    I = Matrix(PauliBitString("I"))
    X = Matrix(PauliBitString("X"))
    Y = Matrix(PauliBitString("Y"))
    Z = Matrix(PauliBitString("Z"))

    σp = (X - 1im*Y) ./ 2
    σm = (X + 1im*Y) ./ 2
    
    N = (I + Z) /2
    M = (I - Z) /2
  
    out = zeros(2^3, 2^3)
    out += kron(N,N,σp) .* sqrt(1)
    out += kron(N,M,σp) .* sqrt(3)
    out += kron(M,N,σp) .* sqrt(5)
    out += kron(M,M,σp) .* sqrt(7)
    out += kron(N,σp,σm) .* sqrt(2)
    out += kron(M,σp,σm) .* sqrt(6)
    out += kron(σp,σm,σm) .* sqrt(4)
    display(out)

    out = zeros(2^4, 2^4)
    out += kron(N,N,N,σp) .* sqrt(1)
    out += kron(N,N,M,σp) .* sqrt(3)
    out += kron(N,M,N,σp) .* sqrt(5)
    out += kron(N,M,M,σp) .* sqrt(7)
    out += kron(M,N,N,σp) .* sqrt(9)
    out += kron(M,N,M,σp) .* sqrt(11)
    out += kron(M,M,N,σp) .* sqrt(13)
    out += kron(M,M,M,σp) .* sqrt(15)
    
    out += kron(N,N,σp,σm) .* sqrt(2)
    out += kron(N,M,σp,σm) .* sqrt(6)
    out += kron(M,N,σp,σm) .* sqrt(10)
    out += kron(M,M,σp,σm) .* sqrt(14)
    
    out += kron(N,σp,σm,σm) .* sqrt(4)
    out += kron(M,σp,σm,σm) .* sqrt(12)
    
    out += kron(σp,σm,σm,σm) .* sqrt(8)
    display(out)

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
            println(prodlist, start + dec_idx*stride)
            push!(rep1, (deepcopy(prodlist), start + dec_idx*stride))
            dec_idx += 1
            
        end
    end

    println(" -------- Initial Representation -------- ")
    # println(rep1)
    for i in rep1
        for j in i[1]
            @printf("%2s ",j)
        end
        # println(i[2])
        @printf("√%i\n",i[2])
    end 


    N = [(.5,"I"), (.5, "Z")]
    M = [(.5,"I"), (-.5, "Z")]
    σp = [(.5,"X"), (-.5im, "Y")]
    σm = [(.5,"X"), (.5im, "Y")]
   

    for i in rep1
    end
    
    return out

    
    
end