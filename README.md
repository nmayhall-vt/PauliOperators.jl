# PauliOperators

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://nmayhall-vt.github.io/PauliOperators.jl/stable/) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://nmayhall-vt.github.io/PauliOperators.jl/dev/)
[![Build Status](https://github.com/nmayhall-vt/PauliOperators.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/nmayhall-vt/PauliOperators.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/nmayhall-vt/PauliOperators.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/nmayhall-vt/PauliOperators.jl)
[![codecov](https://codecov.io/gh/nmayhall-vt/PauliOperators.jl/branch/v2.0/graph/badge.svg?token=8Y2tTRjV4U)](https://codecov.io/gh/nmayhall-vt/PauliOperators.jl)

Simple package providing bitstring representations of tensor products of Pauli operators for efficient manipulations. In this representation, the Pauli string operator is represented as two binary strings, one for x and one for z.

The format used for an arbitrary `Pauli` is as follows: 

$$
\begin{align}
P_n =& i^\theta Z^{z_1} X^{x_1} ⊗ Z^{z_2} X^{x_2} ⊗ ⋯ ⊗ Z^{z_N} X^{x_N}  \\
=& i^\theta \bigotimes_j Z^{z_j} X^{x_j} 
\end{align}
$$ 
where the $z$ and $x$ strings are encoded as the bitwise representation of an integer. Currently, we use `Int128` integers so that have access to relatively large pauli strings, i.e., up to 128 qubits. Using only $X$ and $Z$ operators gives access to $Y$ (or rather $iY$) by simply setting both $x_i$ and $z_i$ to 1.   


## Types
This package provides the following types:

- `FixedPhasePauli`, which contains only the base operator string. Here, `N` is the number of qubits (currently maxed at 128).
```julia
struct FixedPhasePauli{N} <: AbstractPauli{N}
    z::Int128
    x::Int128
end
```
- `Pauli`, which includes both the base operator string, as well as a phase, $\theta$. This allows us to multiply and such, while keep track of the phases.
```julia
struct Pauli{N} <: AbstractPauli{N}
    θ::UInt8
    pauli::FixedPhasePauli{N}
end
```
-  `ScaledPauli`. This allows us to describe a `Pauli` scaled by an arbitrary Float. We may want to parameterize the type here in the future.
```julia
struct ScaledPauli{N} <: AbstractPauli{N}
    coeff::ComplexF64
    pauli::FixedPhasePauli{N}
end
```
-  `PauliSum` which provides a way to define a sum of `ScaledPauli`'s. While there are multiple ways one might choose to do this, here, we use a dictionary, where the `FixedPhasePauli` is the key, and the value is the coefficient. 
```julia
struct PauliSum{N}  
    ops::Dict{FixedPhasePauli{N},ComplexF64}
end
```
-  `ScaledPauliVector` is simply an alias to a vector of `ScaledPauli`'s. This also represents a sum of `ScaledPauli`'s. 
```julia
ScaledPauliVector{N} = Vector{ScaledPauli{N}}
```
- `KetBitString` is a simple bitstring representation of a computational basis state. This is provided to give access to computing things like expectation values and just general wavefunctions. p
```julia
"""
An occupation number vector, up to 128 qubits
"""
struct KetBitString{N} 
    v::Int128
end
```


## Functions
The following functions are overloaded to allow the objects to interact intuitively:
- `*`, `⊗`, `⊕`, `+`

Addition needs a bit more specification, since, unlike `*`, a sum of Pauli's is not a Pauli. As such, we simply choose to return a `PauliSum` type when adding two `AbstractPauli`'s. 
