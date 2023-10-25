module DecomposingPolynomialSystems

using HomotopyContinuation:
    HomotopyContinuation,
    @var,
    System,
    Variable,
    Expression
export @var, System, Variable, Expression
const HC = HomotopyContinuation

import LinearAlgebra
export det, dot, I, nullspace

import GAP

include("utils.jl")
include("monomials.jl")
include("expression_map.jl")
include("sampled_system.jl")
include("scalings.jl")
include("interpolation.jl")
include("deck_transformations.jl")
# include("invariants.jl")
# include("implicitization.jl")
# include("decompose.jl")


end # module
