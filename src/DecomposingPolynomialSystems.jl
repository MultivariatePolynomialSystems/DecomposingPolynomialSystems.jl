module DecomposingPolynomialSystems

using HomotopyContinuation
export @var, System, Variable, Expression

using AbstractAlgebra
using LinearAlgebra
export det, dot

using Combinatorics
using Transducers
using GAP

include("utils/utils.jl")
include("common.jl")
include("symmetries.jl")
include("invariants.jl")
include("implicitization.jl")

struct DecomposablePolynomialSystem
    system::SampledSystem
    is_decomposable::Bool
    deck_transformations::Vector{Vector{Expression}}
    factorizing_maps::Vector{FactorizingMap}
    intermediate_varieties::Vector{SampledSystem}
end

DecomposablePolynomialSystem(system::SampledSystem, is_decomposable::Bool) = DecomposablePolynomialSystem(system, false, [], [], [])
DecomposablePolynomialSystem(system::SampledSystem, deck_transformations::Vector{Vector{Expression}}, factorizing_maps::Vector{FactorizingMap}) = DecomposablePolynomialSystem(system, true, deck_transformations, factorizing_maps, [SampledSystem() for _ in eachindex(factorizing_maps)])

function decompose(F::System, xp0::Tuple{Vector{CC}, Vector{CC}}; degDeck::Int64=1, degFact::Int64=1, degImpl::Int64=1, tol::Float64=1e-5, paramDepDeck::Bool=false, paramDepFact::Bool=false)::DecomposablePolynomialSystem
    F = run_monodromy(F, xp0)
    if length(F.block_partitions) == 0
        return DecomposablePolynomialSystem(F, false)
    end

    n_unknowns = length(xp0[1])
    n_params = length(xp0[2])
    n_sols = size(F.solutions, 2)

    paramDepDeck ? n_vars = n_unknowns + n_params : n_vars = n_unknowns
    n_instances_deck = Int(ceil(2/n_sols*binomial(n_vars + degDeck, n_vars)))

    paramDepFact ? n_vars = n_unknowns + n_params : n_vars = n_unknowns
    n_constraints = min([num_constraints(F.block_partitions[i]) for i in 1:length(F.block_partitions)]...)
    n_instances_fact = Int(ceil(1/n_constraints*binomial(n_vars + degFact, degFact)))

    F = sample_system(F, max(n_instances_deck, n_instances_fact))

    deck_transformations = compute_deck_transformations(F, degree=degDeck, tol=tol, param_dep=paramDepDeck)
    factorizing_maps = compute_factorizing_maps(F, degree=degFact, tol=tol, param_dep=paramDepFact)

    n_new_unknowns = [length(factorizing_maps[i].map) + n_params for i in eachindex(factorizing_maps)]
    n_new_mons = [binomial(n_new_unknowns[i] + degImpl, degImpl) for i in eachindex(factorizing_maps)]
    n_instances_impl = max([Int(ceil(1/length(F.block_partitions[i])*n_new_mons[i])) for i in eachindex(factorizing_maps)]...)
    n_instances_impl = max(n_instances_impl, binomial(n_params + degImpl, degImpl))

    F = sample_system(F, n_instances_impl)

    DPS = DecomposablePolynomialSystem(F, deck_transformations, factorizing_maps)
    for i in eachindex(factorizing_maps)
        factorizing_map = factorizing_maps[i]
        @var y[i, 1:length(factorizing_map.map)]
        DPS.intermediate_varieties[i] = implicitize(F, factorizing_map, F.block_partitions[i], new_vars=y, degree=degImpl)
        DPS.factorizing_maps[i].domain = Ref(DPS.system)
        DPS.factorizing_maps[i].image = Ref(DPS.intermediate_varieties[i])
        DPS.factorizing_maps[i].monodromy_group = action_on_block(F.monodromy_group, F.block_partitions[i])
    end

    return DPS
end

export DecomposablePolynomialSystem
export decompose

end # module
