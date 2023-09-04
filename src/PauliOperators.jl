module PauliOperators

using Printf

include("helpers.jl")
include("type_Pauli.jl")
include("type_ScaledPauli.jl")
include("type_PauliSum.jl")
include("type_KetBitString.jl")
include("operations.jl")
include("transformations.jl")
include("function_mul.jl")
include("function_base.jl")


# Exports
export Pauli
export ScaledPauli
export PauliSum
export KetBitString
export rotate_phase 
export get_phase
export get_coeff
export is_diagonal
export is_hermitian 
export negate 
export commute
export phasefree 
export random_Pauli 
export boson_binary_transformation
export jordan_wigner
export otimes
export clip! 
export ⊗
export ⊕
export nqubits

end

