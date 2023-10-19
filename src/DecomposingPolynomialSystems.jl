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

import Combinatorics
import Transducers
import GAP

include("utils.jl")
include("monomials.jl")
include("expression_map.jl")
include("sampled_system.jl")
include("interpolation.jl")
include("symmetries.jl")
include("invariants.jl")
# include("implicitization.jl")
# include("decompose.jl")


end # module
