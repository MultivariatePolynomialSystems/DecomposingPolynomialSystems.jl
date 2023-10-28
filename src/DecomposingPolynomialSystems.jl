module DecomposingPolynomialSystems

import HomotopyContinuation
const HC = HomotopyContinuation
using HomotopyContinuation.ModelKit
export @var, Variable, Expression, System

using LinearAlgebra: nullspace

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
