using PauliOperators
using LinearAlgebra
using Random
using Test
using BlockDavidson
# using LinearMaps
using BenchmarkTools
using KrylovKit


Random.seed!(1)

N = 3
function build_ops(N)
    H0 = PauliSum(N)
    V = PauliSum(N)
    for z in 0:2^N-1
        H0 += rand(ComplexF64) * FixedPhasePauli{N}(z, 0)
    end
    for z in 0:2^N-1
        for x in 1:2^N-1
            V += rand(ComplexF64) * FixedPhasePauli{N}(z, x)
        end
    end
    H0 = H0 + adjoint(H0)
    V = V + adjoint(V)

    # clip!(H0)
    # clip!(V)
    return H0, V
end
H0, V = build_ops(N)

basis = 0*H0 + 0*V
Asubspace = 0*basis
Bsubspace = 0*basis

Asubspace += H0
Bsubspace += V
Xsubspace = 1im*Bsubspace


function Base.Vector(ps::PauliSum{N}) where {N}
    return Vector{ComplexF64}([i for i in values(ps.ops)])
end

function Base.fill!(ps::PauliSum{N}, v::Vector) where {N}
    length(ps) == length(v) || throw(DimensionMismatch)
    for (i, (pauli,coeff)) in enumerate(ps.ops)
        ps[pauli] = v[i]
    end
end


"""
    mymatvec_subspace(vin)

|i)(i|L|j)(j|x|k) = |i)(i|L|j)x(j) = L(ij)x(j) = tr(Pi' * [A,Pj])x(j)
    = tr(Pi' * A * Pj)x(j) - tr(Pi' * Pj * A)x(j) 
    = sum_k tr(Pi' * Pk * Pj)x(j)a(k) - tr(Pi' * Pj * Pk)x(j)a(k)

if [Pk,Pj]!=0, then {Pk,Pj} = 0, and PkPj = -PjPk

    = tr(Pi' * Pk * Pj)x(j)a(k) + tr(Pi' * Pk * Pj)x(j)a(k)
    = 2*tr(Pi' * Pk * Pj)x(j)a(k) 

    Since A = sum_i ai Pi
    = 2* θ(k*j) x(j) a(k)
"""
function mymatvec_subspace(vin)
    out = deepcopy(basis)
    X = deepcopy(basis)
    fill!(X,vin)
    for (op1, coeff1) in Asubspace.ops #k
        for (op2, coeff2) in X.ops #j
            commute(op1, op2) == false || continue
            prod = op1 * op2
            if haskey(basis, prod)
                out[prod] += 2*coeff1*coeff2
            end
        end
    end
    return Vector(out)
end

"""
    get_full_mat()

|i)(i|L|j)(j| = L(ij) = tr(Pi' * [A,Pj]) = tr(Pi'*A*Pj) - tr(Pi'*Pj*A) = tr(Pj*Pi'*A) - tr(Pi'*Pj*A)

if [Pi',Pj]!=0, then {Pi',Pj} = 0, and Pi'Pj = -PjPi'

    = 2*tr(Pj*Pi'*A)

    Since A = sum_i ai Pi
    = 2*ak θ(ji')

"""
function get_full_mat(basis::PauliSum{N}) where N
    L = zeros(ComplexF64, length(basis), length(basis))

    for (idx1, (op1, coeff1)) in enumerate(basis.ops)   #i
        for (idx2, (op2, coeff2)) in enumerate(basis.ops)   #j
            commute(op1, op2) == false || continue
            prod = op1 * op2
            if haskey(Asubspace, prod)
                L[idx1, idx2] = 2 * Asubspace[prod] * get_phase(prod)
            end

        end
    end
    return L
end

L = get_full_mat(basis)

Aop = LinOpMat{ComplexF64}(mymatvec_subspace, length(basis), false)


btmp1 = Aop * Vector(Xsubspace)
btmp2 = L * Vector(Xsubspace)

all(btmp1 .≈ btmp2) || throw(ErrorException)
Amat = Matrix(Asubspace)
Xmat = Matrix(Xsubspace)
Bmat = Matrix(Bsubspace)

btmp3 = Amat*Xmat - Xmat*Amat 

tmp = deepcopy(basis); fill!(tmp, btmp2); 
norm(Matrix(tmp) - btmp3) < 1e-10 || throw(ErrorException)


alg = GMRES(; maxiter=200, tol=1e-6, verbosity=1)
x, info = linsolve(L, Vector(Bsubspace), ones(ComplexF64, length(Bsubspace)), alg)


fill!(Xsubspace, x)

Amat = Matrix(Asubspace)
Bmat = Matrix(Bsubspace)
Xmat = Matrix(Xsubspace)

error = norm(Amat*Xmat - Xmat*Amat - Bmat)
println(error)