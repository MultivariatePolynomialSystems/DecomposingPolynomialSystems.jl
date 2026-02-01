# Sampling Polynomial Systems

In this Julia package we deal with parametric polynomial systems with finitely many solutions for generic parameters. We use [`HomotopyContinuation.jl`](https://www.juliahomotopycontinuation.org/) to sample such polynomial systems.

## Run monodromy

```@docs
run_monodromy
```

## SampledParametricSystem

`SampledParametricSystem` is a struct type that initially contains a polynomial system, the result of monodromy computations, and the solutions-parameters samples obtained with [`run_monodromy`](@ref) or [`sample!`](@ref).

```@docs
unknowns
parameters
variables
nunknowns
nparameters
nvariables
nsolutions
samples
ninstances
nsamples
monodromy_permutations
block_partitions
aut_permutations
```

## Sample system

```@docs
sample!
```