# DecomposingPolynomialSystems.jl

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://multivariatepolynomialsystems.github.io/DecomposingPolynomialSystems.jl/dev)

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
deck = symmetries_fixing_parameters(F, degree_bound=1, param_dep=false)
```
The object `deck` is a `DeckTransformationGroup` structure that contains 4 deck transformations acting on the unknowns `x₁`, `x₂` of the polynomial system `F`:
```
4-element Vector{Vector{Union{Nothing, Expression}}}:
 [x₁, x₂]
 [(0.0 - 1.0*im)*x₂, (0.0 + 1.0*im)*x₁]
 [(-1.0 + 0.0*im)*x₁, (-1.0 + 0.0*im)*x₂]
 [(0.0 + 1.0*im)*x₂, (0.0 - 1.0*im)*x₁]
```
