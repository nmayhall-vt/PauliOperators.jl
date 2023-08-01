var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = PauliOperators","category":"page"},{"location":"#PauliOperators","page":"Home","title":"PauliOperators","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for PauliOperators.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [PauliOperators]","category":"page"},{"location":"#PauliOperators.KetBitString","page":"Home","title":"PauliOperators.KetBitString","text":"An occupation number vector, up to 128 qubits\n\n\n\n\n\n","category":"type"},{"location":"#PauliOperators.KetBitString-Tuple{Integer, Integer}","page":"Home","title":"PauliOperators.KetBitString","text":"KetBitString(N::Integer, v::Integer)\n\nTBW\n\n\n\n\n\n","category":"method"},{"location":"#PauliOperators.KetBitString-Union{Tuple{Vector{T}}, Tuple{T}} where T<:Integer","page":"Home","title":"PauliOperators.KetBitString","text":"KetBitString(vec::Vector{T}) where T<:Union{Bool, Integer}\n\nTBW\n\n\n\n\n\n","category":"method"},{"location":"#PauliOperators.PauliBitString","page":"Home","title":"PauliOperators.PauliBitString","text":"In this representation, the Pauli string operator is represented as two binary strings, one for x and one for z.\n\nThe format is as follows: \n\n(i)^\theta Z^z1 X^x1 ⊗ Z^z2 X^x2 ⊗ ⋯ ⊗ Z^zN X^xN\n\nProducts of operators simply concatonate the left and right strings separately. For example, \n\nXYZIY = 11001|01101\n\nTo create a Y operator, bits in the same locations in z and x should be on.  This means that we have a phase to keep track of because Z^1*X^1 = iY.  As such, we end up working with operators of the form:\n\n(i)^\theta σ1 ⊗ σ2 ⊗ σ3 ⊗ ⋯ ⊗ σN,\n\nwhere,\n\nσ ∈ {X, iY, Z, I}\n\n\n\n\n\n","category":"type"},{"location":"#PauliOperators.PauliBitString-Tuple{Integer}","page":"Home","title":"PauliOperators.PauliBitString","text":"PauliBitString(N::Integer; X=[], Y=[], Z=[])\n\nconstructor for creating PauliBoolVec by specifying the qubits where each X, Y, and Z gates exist \n\n\n\n\n\n","category":"method"},{"location":"#PauliOperators.PauliBitString-Tuple{String}","page":"Home","title":"PauliOperators.PauliBitString","text":"PauliBitString(str::String)\n\nTBW\n\n\n\n\n\n","category":"method"},{"location":"#PauliOperators.PauliBitString-Union{Tuple{I}, Tuple{I, I}} where I<:Integer","page":"Home","title":"PauliOperators.PauliBitString","text":"PauliBitString(z::I, x::I) where I<:Integer\n\nTBW\n\n\n\n\n\n","category":"method"},{"location":"#Base.:*-Union{Tuple{N}, Tuple{PauliBitString{N}, PauliBitString{N}}} where N","page":"Home","title":"Base.:*","text":"Base.:*(p1::PauliBitString{N}, p2::PauliBitString{N}) where {N}\n\nMultiply two PauliBitString's together\n\n\n\n\n\n","category":"method"},{"location":"#Base.:==-Union{Tuple{N}, Tuple{PauliBitString{N}, PauliBitString{N}}} where N","page":"Home","title":"Base.:==","text":"Base.:(==)(p1::PauliBitString{N}, p2::PauliBitString{N}) where {N}\n\nCheck if they are equal, return true or false\n\n\n\n\n\n","category":"method"},{"location":"#Base.show-Union{Tuple{N}, Tuple{IO, KetBitString{N}}} where N","page":"Home","title":"Base.show","text":"Base.show(io::IO, P::PauliBitString{N}) where N\n\nTBW\n\n\n\n\n\n","category":"method"},{"location":"#Base.show-Union{Tuple{N}, Tuple{IO, PauliBitString{N}}} where N","page":"Home","title":"Base.show","text":"Base.show(io::IO, P::PauliMask)\n\nTBW\n\n\n\n\n\n","category":"method"},{"location":"#Base.string-Union{Tuple{KetBitString{N}}, Tuple{N}} where N","page":"Home","title":"Base.string","text":"Base.string(p::KetBitString{N}) where N\n\nDisplay, y = iY\n\n\n\n\n\n","category":"method"},{"location":"#Base.string-Union{Tuple{PauliBitString{N}}, Tuple{N}} where N","page":"Home","title":"Base.string","text":"Base.display(p::PauliBitString)\n\nDisplay, y = iY\n\n\n\n\n\n","category":"method"},{"location":"#PauliOperators.commute-Union{Tuple{N}, Tuple{PauliBitString{N}, PauliBitString{N}}} where N","page":"Home","title":"PauliOperators.commute","text":"commute(p1::PauliBitString{N}, p2::PauliBitString{N}) where {N}\n\nCheck if they commute, return true or false\n\n\n\n\n\n","category":"method"},{"location":"#PauliOperators.expectation_value_sign-Union{Tuple{N}, Tuple{PauliBitString{N}, KetBitString{N}}} where N","page":"Home","title":"PauliOperators.expectation_value_sign","text":"expectation_value_sign(p::PauliBitString{N}, ket::Vector{Bool}) where N\n\ncompute expectation value of PauliBitString o for a product state ket\n\n\n\n\n\n","category":"method"},{"location":"#PauliOperators.get_phase-Tuple{PauliBitString}","page":"Home","title":"PauliOperators.get_phase","text":"get_phase(p::PauliBitString)\n\nReturn the phase of the PauliBitString, i^θ\n\n\n\n\n\n","category":"method"},{"location":"#PauliOperators.is_diagonal-Tuple{PauliBitString}","page":"Home","title":"PauliOperators.is_diagonal","text":"is_diagonal(p::PauliBitString)\n\nCheck if operator is diagonal in the computational (z) basis. E.g., does this operator consist of only I and/or Z?\n\n\n\n\n\n","category":"method"},{"location":"#PauliOperators.negate-Union{Tuple{PauliBitString{N}}, Tuple{N}} where N","page":"Home","title":"PauliOperators.negate","text":"negate(p::PauliBitString)\n\nMultiply p by -1\n\n\n\n\n\n","category":"method"},{"location":"#PauliOperators.random_PauliBitString-Tuple{Any}","page":"Home","title":"PauliOperators.random_PauliBitString","text":"random_PauliBitString(N)\n\nTBW\n\n\n\n\n\n","category":"method"},{"location":"#PauliOperators.rotate_phase-Union{Tuple{N}, Tuple{PauliBitString{N}, Integer}} where N","page":"Home","title":"PauliOperators.rotate_phase","text":"rotate_phase(p::PauliBitString{N}, θ::Integer) where N\n\nRotate phase in units of π/2. In otherwords, multiply the phase by i^θ. E.g., mutliplication by -1 is obtained with θ=2.\n\n\n\n\n\n","category":"method"}]
}