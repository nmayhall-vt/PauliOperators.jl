module PauliOperators

    using Printf
    using LinearAlgebra
   

    include("helpers.jl")
    include("type_PauliBasis.jl")
    include("type_Pauli.jl")
    include("type_PauliSum.jl")
    include("type_Ket.jl")
    include("type_KetSum.jl")
    include("type_DyadBasis.jl")
    include("type_Dyad.jl")
    include("type_DyadSum.jl")
    include("multiplication.jl")
    include("addition.jl")
    include("conversions.jl")

    const ⊗ = otimes
    const ⊕ = osum

    export Pauli
    export PauliBasis
    export PauliSum
    export Ket
    export Bra
    export DyadBasis
    export Dyad
    export DyadSum
    export KetSum
    export clip! 
    export ⊗
    export ⊕
    export expectation_value

    export symplectic_phase
    export coeff
end
