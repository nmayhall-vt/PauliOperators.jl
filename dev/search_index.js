var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = PauliOperators","category":"page"},{"location":"#PauliOperators","page":"Home","title":"PauliOperators","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for PauliOperators.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [PauliOperators]","category":"page"},{"location":"#Base.Matrix-Union{Tuple{Pauli{N}}, Tuple{N}} where N","page":"Home","title":"Base.Matrix","text":"Base.Matrix(p::Pauli{N}) where N\n\nCreate dense matrix representation \n\n\n\n\n\n","category":"method"},{"location":"#PauliOperators.KetBitString","page":"Home","title":"PauliOperators.KetBitString","text":"An occupation number vector, up to 128 qubits\n\n\n\n\n\n","category":"type"},{"location":"#PauliOperators.KetBitString-Tuple{Integer, Integer}","page":"Home","title":"PauliOperators.KetBitString","text":"KetBitString(N::Integer, v::Integer)\n\nTBW\n\n\n\n\n\n","category":"method"},{"location":"#PauliOperators.KetBitString-Union{Tuple{Vector{T}}, Tuple{T}} where T<:Integer","page":"Home","title":"PauliOperators.KetBitString","text":"KetBitString(vec::Vector{T}) where T<:Union{Bool, Integer}\n\nTBW\n\n\n\n\n\n","category":"method"},{"location":"#PauliOperators.Pauli","page":"Home","title":"PauliOperators.Pauli","text":"In this representation, the Pauli string operator is represented as two binary strings, one for x and one for z.\n\nThe format is as follows: \n\n(i)^\theta Z^z1 X^x1 ⊗ Z^z2 X^x2 ⊗ ⋯ ⊗ Z^zN X^xN\n\nProducts of operators simply concatonate the left and right strings separately. For example, \n\nXYZIY = 11001|01101\n\nTo create a Y operator, bits in the same locations in z and x should be on.  This means that we have a phase to keep track of because Z^1*X^1 = iY.  As such, we end up working with operators of the form:\n\n(i)^\theta σ1 ⊗ σ2 ⊗ σ3 ⊗ ⋯ ⊗ σN,\n\nwhere,\n\nσ ∈ {X, iY, Z, I}\n\n\n\n\n\n","category":"type"},{"location":"#PauliOperators.Pauli-Tuple{Integer}","page":"Home","title":"PauliOperators.Pauli","text":"Pauli(N::Integer; X=[], Y=[], Z=[])\n\nconstructor for creating PauliBoolVec by specifying the qubits where each X, Y, and Z gates exist \n\n\n\n\n\n","category":"method"},{"location":"#PauliOperators.Pauli-Tuple{String}","page":"Home","title":"PauliOperators.Pauli","text":"Pauli(str::String)\n\nTBW\n\n\n\n\n\n","category":"method"},{"location":"#PauliOperators.Pauli-Union{Tuple{I}, Tuple{I, I, Any}} where I<:Integer","page":"Home","title":"PauliOperators.Pauli","text":"Pauli(z::I, x::I) where I<:Integer\n\nTBW\n\n\n\n\n\n","category":"method"},{"location":"#PauliOperators.PauliSum","page":"Home","title":"PauliOperators.PauliSum","text":"ops::Dict{Pauli{N},ComplexF64}\n\nA collection of Paulis, joined by addition\n\n\n\n\n\n","category":"type"},{"location":"#PauliOperators.PauliSum-Tuple{Any}","page":"Home","title":"PauliOperators.PauliSum","text":"PauliSum(N)\n\nTBW\n\n\n\n\n\n","category":"method"},{"location":"#Base.:*-Union{Tuple{N}, Tuple{PauliSum{N}, Number}} where N","page":"Home","title":"Base.:*","text":"Base.:*(ps::PauliSum{N}, a::Number) where {N}\n\nTBW\n\n\n\n\n\n","category":"method"},{"location":"#Base.:*-Union{Tuple{N}, Tuple{PauliSum{N}, PauliSum{N}}} where N","page":"Home","title":"Base.:*","text":"Base.-(p1::PauliSum{N}, p2::PauliSum{N}) where {N}\n\nMultiply two PauliSums. \n\n\n\n\n\n","category":"method"},{"location":"#Base.:*-Union{Tuple{N}, Tuple{Pauli{N}, Pauli{N}}} where N","page":"Home","title":"Base.:*","text":"Base.:*(p1::Pauli{N}, p2::Pauli{N}) where {N}\n\nMultiply two Pauli's together\n\n\n\n\n\n","category":"method"},{"location":"#Base.:+-Union{Tuple{N}, Tuple{PauliSum{N}, PauliSum{N}}} where N","page":"Home","title":"Base.:+","text":"Base.+(p1::PauliSum{N}, p2::PauliSum{N}) where {N}\n\nAdd two PauliSums. \n\n\n\n\n\n","category":"method"},{"location":"#Base.:+-Union{Tuple{N}, Tuple{PauliSum{N}, Pauli{N}}} where N","page":"Home","title":"Base.:+","text":"Base.:+(ps::PauliSum{N}, p::Pauli{N}) where {N}\n\nAdd a Pauli to a PauliSum. \n\n\n\n\n\n","category":"method"},{"location":"#Base.:+-Union{Tuple{N}, Tuple{Pauli{N}, Pauli{N}}} where N","page":"Home","title":"Base.:+","text":"Base.:+(p1::Pauli{N}, p2::Pauli{N}) where {N}\n\nAdd two Pauli's together. This returns a PauliSum\n\n\n\n\n\n","category":"method"},{"location":"#Base.:--Union{Tuple{N}, Tuple{PauliSum{N}, PauliSum{N}}} where N","page":"Home","title":"Base.:-","text":"Base.-(p1::PauliSum{N}, p2::PauliSum{N}) where {N}\n\nSubtract two PauliSums. \n\n\n\n\n\n","category":"method"},{"location":"#Base.:--Union{Tuple{N}, Tuple{PauliSum{N}, Pauli{N}}} where N","page":"Home","title":"Base.:-","text":"Base.:-(ps::PauliSum{N}, p::Pauli{N}) where {N}\n\nAdd a Pauli to a PauliSum. \n\n\n\n\n\n","category":"method"},{"location":"#Base.:--Union{Tuple{N}, Tuple{Pauli{N}, Pauli{N}}} where N","page":"Home","title":"Base.:-","text":"Base.:-(p1::Pauli{N}, p2::Pauli{N}) where {N}\n\nSubtract two Pauli's. This returns a PauliSum\n\n\n\n\n\n","category":"method"},{"location":"#Base.:<-Union{Tuple{N}, Tuple{Pauli{N}, Pauli{N}}} where N","page":"Home","title":"Base.:<","text":"Base.:(<)(p1::Pauli{N}, p2::Pauli{N}) where {N}\n\nCheck if p1 < p2\n\n\n\n\n\n","category":"method"},{"location":"#Base.:==-Union{Tuple{N}, Tuple{Pauli{N}, Pauli{N}}} where N","page":"Home","title":"Base.:==","text":"Base.:(==)(p1::Pauli{N}, p2::Pauli{N}) where {N}\n\nCheck if they are equal, return true or false\n\n\n\n\n\n","category":"method"},{"location":"#Base.:>-Union{Tuple{N}, Tuple{Pauli{N}, Pauli{N}}} where N","page":"Home","title":"Base.:>","text":"Base.:(>)(p1::Pauli{N}, p2::Pauli{N}) where {N}\n\nCheck if p1 > p2\n\n\n\n\n\n","category":"method"},{"location":"#Base.Multimedia.display-Tuple{PauliSum}","page":"Home","title":"Base.Multimedia.display","text":"Base.display(ps::PauliSum)\n\nTBW\n\n\n\n\n\n","category":"method"},{"location":"#Base.adjoint-Union{Tuple{PauliSum{N}}, Tuple{N}} where N","page":"Home","title":"Base.adjoint","text":"LinearAlgebra.adjoint(ps::PauliSum{N}) where N\n\nTBW\n\n\n\n\n\n","category":"method"},{"location":"#Base.hash-Union{Tuple{N}, Tuple{Pauli{N}, UInt64}} where N","page":"Home","title":"Base.hash","text":"Base.hash(p::Pauli{N}, h::UInt) where N\n\nCreate a hash for a Pauli. Because we want to collect matching operators,      with different phases, we don't actually put the phase in the hash\n\n\n\n\n\n","category":"method"},{"location":"#Base.show-Union{Tuple{N}, Tuple{IO, KetBitString{N}}} where N","page":"Home","title":"Base.show","text":"Base.show(io::IO, P::Pauli{N}) where N\n\nTBW\n\n\n\n\n\n","category":"method"},{"location":"#Base.show-Union{Tuple{N}, Tuple{IO, Pauli{N}}} where N","page":"Home","title":"Base.show","text":"Base.show(io::IO, P::PauliMask)\n\nTBW\n\n\n\n\n\n","category":"method"},{"location":"#Base.string-Union{Tuple{KetBitString{N}}, Tuple{N}} where N","page":"Home","title":"Base.string","text":"Base.string(p::KetBitString{N}) where N\n\nDisplay, y = iY\n\n\n\n\n\n","category":"method"},{"location":"#Base.string-Union{Tuple{Pauli{N}}, Tuple{N}} where N","page":"Home","title":"Base.string","text":"Base.display(p::Pauli)\n\nDisplay, y = iY\n\n\n\n\n\n","category":"method"},{"location":"#Base.sum!-Union{Tuple{N}, Tuple{PauliSum{N}, PauliSum{N}}} where N","page":"Home","title":"Base.sum!","text":"Base.sum!(p1::PauliSum{N}, p2::PauliSum{N}) where {N}\n\nAdd two PauliSums. \n\n\n\n\n\n","category":"method"},{"location":"#Base.sum!-Union{Tuple{N}, Tuple{PauliSum{N}, Pauli{N}}} where N","page":"Home","title":"Base.sum!","text":"Base.sum!(p1::PauliSum{N}, p2::Pauli{N}) where {N}\n\nAdd a Pauli to a PauliSum. \n\n\n\n\n\n","category":"method"},{"location":"#LinearAlgebra.adjoint!-Union{Tuple{PauliSum{N}}, Tuple{N}} where N","page":"Home","title":"LinearAlgebra.adjoint!","text":"LinearAlgebra.adjoint!(ps::PauliSum{N}) where N\n\nTBW\n\n\n\n\n\n","category":"method"},{"location":"#PauliOperators.boson_binary_transformation-Tuple{Any}","page":"Home","title":"PauliOperators.boson_binary_transformation","text":"boson_binary_transformation(nqubits)\n\nTBW\n\n\n\n\n\n","category":"method"},{"location":"#PauliOperators.commute-Union{Tuple{N}, Tuple{Pauli{N}, Pauli{N}}} where N","page":"Home","title":"PauliOperators.commute","text":"commute(p1::Pauli{N}, p2::Pauli{N}) where {N}\n\nCheck if they commute, return true or false\n\n\n\n\n\n","category":"method"},{"location":"#PauliOperators.expectation_value_sign-Union{Tuple{N}, Tuple{Pauli{N}, KetBitString{N}}} where N","page":"Home","title":"PauliOperators.expectation_value_sign","text":"expectation_value_sign(p::Pauli{N}, ket::Vector{Bool}) where N\n\ncompute expectation value of Pauli o for a product state ket\n\n\n\n\n\n","category":"method"},{"location":"#PauliOperators.get_phase-Tuple{Pauli}","page":"Home","title":"PauliOperators.get_phase","text":"get_phase(p::Pauli)\n\nReturn the phase of the Pauli, i^θ\n\n\n\n\n\n","category":"method"},{"location":"#PauliOperators.is_diagonal-Tuple{Pauli}","page":"Home","title":"PauliOperators.is_diagonal","text":"is_diagonal(p::Pauli)\n\nCheck if operator is diagonal in the computational (z) basis. E.g., does this operator consist of only I and/or Z?\n\n\n\n\n\n","category":"method"},{"location":"#PauliOperators.mult!-Tuple{PauliSum, Number}","page":"Home","title":"PauliOperators.mult!","text":"Base.mult!(ps::PauliSum{N}, a::Number)\n\nTBW\n\n\n\n\n\n","category":"method"},{"location":"#PauliOperators.negate-Union{Tuple{Pauli{N}}, Tuple{N}} where N","page":"Home","title":"PauliOperators.negate","text":"negate(p::Pauli)\n\nMultiply p by -1\n\n\n\n\n\n","category":"method"},{"location":"#PauliOperators.random_Pauli-Tuple{Any}","page":"Home","title":"PauliOperators.random_Pauli","text":"random_Pauli(N)\n\nTBW\n\n\n\n\n\n","category":"method"},{"location":"#PauliOperators.rotate_phase-Union{Tuple{N}, Tuple{Pauli{N}, Integer}} where N","page":"Home","title":"PauliOperators.rotate_phase","text":"rotate_phase(p::Pauli{N}, θ::Integer) where N\n\nRotate phase in units of π/2. In otherwords, multiply the phase by i^θ. E.g., mutliplication by -1 is obtained with θ=2.\n\n\n\n\n\n","category":"method"},{"location":"#PauliOperators.subtract!-Union{Tuple{N}, Tuple{PauliSum{N}, Pauli{N}}} where N","page":"Home","title":"PauliOperators.subtract!","text":"Base.sum!(p1::PauliSum{N}, p2::Pauli{N}) where {N}\n\nAdd a Pauli to a PauliSum. \n\n\n\n\n\n","category":"method"}]
}
