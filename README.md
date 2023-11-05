# DecomposingPolynomialSystems.jl

[![License](https://img.shields.io/badge/License-MIT-yellow?style=flat&label=License&color=%23FFA500%09)](https://github.com/MultivariatePolynomialSystems/DecomposingPolynomialSystems.jl/blob/main/LICENSE)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://multivariatepolynomialsystems.github.io/DecomposingPolynomialSystems.jl/dev)

DecomposingPolynomialSystems.jl is a Julia package that computes the symmetries that fix the parameters (specifically, the group of deck transformations) of a parametric polynomial system with finitely many solutions for generic parameters with a view towards decomposing the given polynomial system.

## Installation

Enter the Pkg REPL by pressing `]` from the Julia REPL and then type
```julia
add https://github.com/MultivariatePolynomialSystems/DecomposingPolynomialSystems.jl.git
```
To get back to the Julia REPL, press backspace.

## Usage
### Computing symmetries
```julia
using DecomposingPolynomialSystems
@var x[1:2] p[1:2]
F = System([x[1]^2 - x[2]^2 - p[1], 2*x[1]*x[2] - p[2]]; variables=x, parameters=p)
symmetries_fixing_parameters(F; degree_bound=1, param_dep=false)
```
The result of the last command is the object of type `DeckTransformationGroup` that contains 4 deck transformations acting on the unknowns `x₁`, `x₂` of the polynomial system `F`:
```
DeckTransformationGroup of order 4
 structure: C2 x C2
 action:
  1st map:
   x₁ ↦ x₁
   x₂ ↦ x₂
  2nd map:
   x₁ ↦ (0.0 + 1.0*im)*x₂
   x₂ ↦ (0.0 - 1.0*im)*x₁
  3rd map:
   x₁ ↦ (0.0 - 1.0*im)*x₂
   x₂ ↦ (0.0 + 1.0*im)*x₁
  4th map:
   x₁ ↦ (-1.0 + 0.0*im)*x₁
   x₂ ↦ (-1.0 + 0.0*im)*x₂
```
### Computing invariants
TBA...
### Decomposition
TBA...