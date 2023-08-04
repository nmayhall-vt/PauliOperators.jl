"""
In this representation, the Pauli string operator is represented as two binary strings, one for x and one for z.
However, unlike `Pauli` types, the `PauliPhaseless` type does not have a phase to keep track of. 

The format is as follows: 
    
    Z^z1 X^x1 ⊗ Z^z2 X^x2 ⊗ ⋯ ⊗ Z^zN X^xN  
    
Products of operators simply concatonate the left and right strings separately. For example, 

    XYZIY = 11001|01101


To create a Y operator, bits in the same locations in `z` and `x` should be on. 
This means that our multiplication is not complete, as we should have a phase to keep track of because Z^1*X^1 = iY.
However, for `PauliPhaseless`, we have Z*X=Y. 
"""
struct Pauliless{N} <: Integer
    z::Int128
    x::Int128
end