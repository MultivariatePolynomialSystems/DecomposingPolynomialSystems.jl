# SampledSystem

`SampledSystem` is a struct type that initially contains a polynomial system, the result of monodromy computations, and the solutions-parameters samples obtained by monodromy or further by additional tracking. This is the return type of [`run_monodromy`](@ref) and [`sample!`](@ref).

## Run monodromy

```@docs
run_monodromy
```

## Sample system

```@docs
sample!
```