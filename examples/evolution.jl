using PauliOperators
using Printf

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
    value_clip!(O, thresh=ϵ)
    return 
end

function evolve_matrix!(O::Matrix, G::Matrix, dt) where {N}
    U = exp(-G * dt * 1im / 2)
    O .= U' * O * U
    return 
end

function run_evolve(N)
    H = heisenberg_1D(N, 1, 1, 1, x=.01, y=.02, z=.03)

    O = PauliSum(N)
    O += Pauli(N, Z=[1])
    Omat = Matrix(O)
    psi = Ket{N}(0)
    print(psi)
    
    print(" Hamiltonian:\n")
    display(H)

    print(" Observable:\n")
    display(O)

    generators = [i*j for (i,j) in H ]

    dt = .1

    for i in 1:1000

        for Gi in generators
            evolve_clip!(O, Gi, dt, ϵ=1e-5)
            evolve_matrix!(Omat, Matrix(Gi), dt)
        end

        val = expectation_value(O, psi)
        ref = tr(Omat * Matrix(DyadBasis{N}(psi,psi')))

        @printf(" val = %12.8f ref = %12.8f err = %12.8f len: %4i\n", val, real(ref), real(val-ref), length(O))
    end
end

run_evolve(4)