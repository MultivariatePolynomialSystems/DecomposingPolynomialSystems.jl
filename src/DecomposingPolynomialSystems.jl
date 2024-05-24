module DecomposingPolynomialSystems

import HomotopyContinuation
const HC = HomotopyContinuation
using HomotopyContinuation.ModelKit
export @var, Variable, Expression, System

using SparseArrays: SparseVector, SparseMatrixCSC, spzeros, AbstractSparseVector, findnz, sparse
using Combinatorics: partitions, multiset_permutations, combinations
using LinearAlgebra: nullspace
using UnPack: @unpack

include("utils.jl")
include("sampled_system.jl")
include("monomials.jl")
# include("expression_map.jl")
include("scalings.jl")
include("interpolation.jl")
include("deck_transformations.jl")
# include("invariants.jl")
# include("implicitization.jl")
# include("decompose.jl")


end # module
