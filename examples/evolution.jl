using PauliOperators
using Printf
using LinearAlgebra
using Random

function heisenberg_1D(N, Jx, Jy, Jz; x=0, y=0, z=0)
    H = PauliSum(N, Float64)
    for i in 0:N-1
        H += -2*Jx * Pauli(N, X=[i+1,(i+1)%(N)+1])
        H += -2*Jy * Pauli(N, Y=[i+1,(i+1)%(N)+1])
        H += -2*Jz * Pauli(N, Z=[i+1,(i+1)%(N)+1])
    end 
    for i in 1:N
        H += x * Pauli(N, X=[i])
        H += y * Pauli(N, Y=[i])
        H += z * Pauli(N, Z=[i])
    end 
    return H
end 

function value_clip!(ps::PauliSum{N}; thresh=1e-16) where {N}
    filter!(p->abs(p.second) > thresh, ps)
    # filter!(p->abs(p.second)/weight(p.first) > thresh, ps.ops)
end

weight(p) = count_ones(p.x | p.z)
# function dampen!(ps::PauliSum{N}, n, alpha) where {N}
#     for (pi,ci) in ps
#         if weight(pi) > n
#             ps[pi] = ci * exp(- alpha * (weight(pi)-n))
#         end
#     end
#     # map!(p->p.second/weight(p.first)^alpha , ps)
# end

function weight_clip!(ps::PauliSum{N}, n) where {N}
    filter!(p->weight(p.first) <= n, ps)
end

function LinearAlgebra.norm(ps::PauliSum)
    l2 = 0
    for (pi,ci) in ps
        l2 += ci*ci'
    end
    return sqrt(l2)
end

function evolve!(O::PauliSum{N,T}, G::Pauli{N}, dt) where {N,T}
    _cos = cos(dt*coeff(G))
    _sin = 1im*sin(dt*coeff(G))
    
    sin_branch = PauliSum(N)
    for (Oi, ci) in O
        if PauliOperators.commute(Oi,PauliBasis(G)) == false
            O[Oi] *= _cos
            sum!(sin_branch, ci * _sin * PauliBasis(G) * Oi)
        end
    end
    sum!(O, sin_branch)
    return 
end

function evolve(O, G, dt)
    O2 = deepcopy(O)
    evolve!(O2, G, dt)
    return O2
end

function evolve_clip!(O::PauliSum{N,T}, G::Pauli{N}, dt; ϵ=1e-3) where {N,T}
    _cos = cos(dt*coeff(G))
    _sin = -1im*sin(dt*coeff(G))
    # trig_min = minimum(_cos, _sin)
    
    sin_branch = PauliSum(N)
    for (Oi, ci) in O
        if PauliOperators.commute(Oi,PauliBasis(G)) == false
            O[Oi] *= _cos
            sum!(sin_branch, ci * _sin * Oi * PauliBasis(G))
        end
    end
    sum!(O, sin_branch)
    weight_clip!(O, 4)
    value_clip!(O, thresh=ϵ)
    return 
end

function evolve_over_rotate!(O::PauliSum{N,T}, G::Pauli{N}, dt; ϵ=1e-3) where {N,T}
    _cos = cos(dt*coeff(G))
    _sin = sin(dt*coeff(G))
    trig_min = min(abs(_cos), abs(_sin))
    
    sin_branch = PauliSum(N)
    for (Oi, ci) in O
        if PauliOperators.commute(Oi,PauliBasis(G)) == false
            if abs(ci) < ϵ
               continue 
            # if abs(trig_min * ci) < ϵ
                # continue
                # transformation will result in a clip
                if abs(_sin) > abs(_cos)
                    # ci -> pi/2
                    O[Oi] = 0
                    sum!(sin_branch, 1im * ci * PauliBasis(G) * Oi)
                else
                    # ci -> 0 
                end
            else
                O[Oi] *= _cos
                sum!(sin_branch, 1im * ci * _sin * PauliBasis(G) * Oi)
            end
        end
    end
    sum!(O, sin_branch)
    dampen!(O, 0)
    # value_clip!(O, thresh=ϵ)
    return 
end

function evolve_matrix!(O::Matrix, G::Matrix, dt)
    U = exp(-G * dt * 1im / 2)
    O .= U' * O * U
    return 
end

function frame_rotation(O, G, rho::DyadSum, V::PauliBasis{N}) where N
    Onew = evolve(O, Pauli(V), π/2)
    Gnew = []
    for Gi in G
        if commute(PauliBasis(Gi), V)
            push!(Gnew, Gi)
        else
            push!(Gnew, 1im * V * Gi)
        end
    end
    
    # Gnew = [evolve(PauliSum(Gi), Pauli(V), π/2) for Gi in G] 
    # Hnew = evolve(H, V*π/2)
    # Onew = 1im * O * V
    # Hnew = 1im * H * V
    # Gnew = [1im * gi * V for gi in generators]
    UV = (Pauli{N}(1,0,0) + 1im*V) * (1/sqrt(2))
    display(UV)
    rhonew = UV * rho * UV'
    # rhonew = (rho + V * rho * V - 1im * rho * V + 1im * V * rho ) * .5 
    # println(tr(Matrix(rho)))
    println(" New trace: ", tr(Matrix(rhonew)))
    # display(Matrix(rhonew))
    # print(sgn)
    return Onew, Gnew, rhonew
    # return Onew, Hnew, rhonew
end

function run_evolve(N)
    H = heisenberg_1D(N, 1, 1, 1, x=.01, y=.02, z=.03)
    

    O = PauliSum(N)
    O += Pauli(N, Z=[1])
    rho = DyadSum(Ket{N}(0), T=ComplexF64)
    display(rho)

    generators = [i*j for (i,j) in H ]
   
    # Rotate Frame
    Random.seed!(2)
    V = rand(PauliBasis{N}) #π/2
    O, generators, rho = frame_rotation(O, generators, rho, V)
 
    println("generators:")
    for i in generators
        display(i)
    end
    # @show typeof(generators)
    # return 

    # return
    println(" New ψ: ")
    display(rho)

    println("random V: ")
    display(V)

    Omat = Matrix(O)
    ρmat = Matrix(rho)
    
    print(" Hamiltonian:\n")
    display(H)

    print(" Observable:\n")
    display(O)
    


    
    dt = .5

    for i in 1:10

        for Gi in generators
            evolve_clip!(O, Gi, dt, ϵ=1e-3)
            # evolve_over_rotate!(O, Gi, dt, ϵ=1e-2)
            evolve_matrix!(Omat, Matrix(Gi), dt)
        end

        val = expectation_value(O, rho)
        ref = tr(Omat * ρmat)

        # println("ref: ", ref)

        @printf("%-4i val = %12.8f ref = %12.8f err = %12.8f len: %4i l2: %12.8f\n", i, real(val), real(ref), real(val-ref), length(O), norm(O))
    end

    # display(O)
end

run_evolve(6)