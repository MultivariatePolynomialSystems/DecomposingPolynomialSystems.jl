module DecomposingPolynomialSystems

import HomotopyContinuation
const HC = HomotopyContinuation
using HomotopyContinuation.ModelKit
export @var, Variable, Expression

using SparseArrays: SparseVector, SparseMatrixCSC, spzeros, AbstractSparseVector, findnz, sparse
using Combinatorics: partitions, multiset_permutations, combinations
using LinearAlgebra: nullspace, I, det
export det
using UnPack: @unpack

include("utils.jl")
include("param_system.jl")
include("sampling/sampled_system.jl")
include("monomials.jl")
include("scalings.jl")
include("interpolation.jl")
include("deck_transformations.jl")

end # module
