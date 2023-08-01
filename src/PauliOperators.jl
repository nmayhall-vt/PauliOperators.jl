module PauliOperators

using Printf




include("helpers.jl")
include("type_PauliBitString.jl")
include("type_KetBitString.jl")
include("operations.jl")


# Exports
export PauliBitString
export KetBitString
export rotate_phase 
export get_phase
export is_diagonal
export negate 

end
