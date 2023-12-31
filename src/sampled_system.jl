export SampledSystem,
    MonodromyInfo,
    Samples,
    run_monodromy,
    sample!,
    unknowns,
    parameters,
    variables,
    nunknowns,
    nparameters,
    nvariables,
    nsolutions,
    nsamples,
    ninstances,
    samples,
    monodromy_permutations,
    block_partitions,
    deck_permutations

using HomotopyContinuation: Result, MonodromyResult, nsolutions, ntracked, is_success, solution
using HomotopyContinuation: ParameterHomotopy, Tracker, track

const MONODROMY_SOLVE_REF = "https://www.juliahomotopycontinuation.org/HomotopyContinuation.jl/stable/monodromy/"
const SOLVE_REF = "https://www.juliahomotopycontinuation.org/HomotopyContinuation.jl/stable/solve/"

struct MonodromyInfo
    n_solutions::Int
    monodromy_permutations::Vector{Vector{Int}}
    block_partitions::Vector{Vector{Vector{Int}}}
    deck_permutations::Vector{Vector{Int}}
end

MonodromyInfo() = MonodromyInfo(1, [], [], [])

struct Samples
    solutions::Array{ComplexF64, 3}  # n_unknowns x n_sols x n_instances
    parameters::Array{ComplexF64, 2}  # n_params x n_instances
end

Samples(
    sols::Matrix{ComplexF64},
    params::Vector{ComplexF64}
) = Samples(reshape(sols, size(sols)..., 1), reshape(params, :, 1))

HC.nsolutions(samples::Samples) = size(samples.solutions, 2)
ninstances(samples::Samples) = size(samples.parameters, 2)
nsamples(samples::Samples) = nsolutions(samples)*ninstances(samples)

mutable struct SampledSystem
    system::System
    mon_info::MonodromyInfo
    samples::Dict{Vector{Int}, Samples} # key: ids of solution paths
end

unknowns(F::System) = HC.variables(F)

"""
    unknowns(F::SampledSystem) -> Vector{Variable}

Returns the vector of unknowns of `F`.
"""
unknowns(F::SampledSystem) = unknowns(F.system)

"""
    parameters(F::SampledSystem) -> Vector{Variable}

Returns the vector of parameters of `F`.
"""
HC.parameters(F::SampledSystem) = parameters(F.system)

variables(F::System) = vcat(unknowns(F), parameters(F))  # does a different thing than HC.variables

"""
    variables(F::SampledSystem) -> Vector{Variable}

Returns the concatenated vector of unknowns and parameters of `F`.
"""
variables(F::SampledSystem) = variables(F.system)

"""
    nunknowns(F::SampledSystem) -> Int

Returns the number of unknowns of `F`.
"""
nunknowns(F::SampledSystem) = length(unknowns(F))  # equivalent to HC.nvariables

"""
    nparameters(F::SampledSystem) -> Int

Returns the number of parameters of `F`.
"""
HC.nparameters(F::SampledSystem) = length(parameters(F))

"""
    nvariables(F::SampledSystem) -> Int

Returns the number of variables of `F`.
"""
nvariables(F::SampledSystem) = length(variables(F))  # doesn't extend HC.nvariables, does a different thing

"""
    nsolutions(F::SampledSystem) -> Int

Returns the number of solutions of `F` obtained by [`run_monodromy`](@ref) method.
"""
HC.nsolutions(F::SampledSystem) = F.mon_info.n_solutions

"""
    samples(F::SampledSystem) -> Dict{Vector{Int}, Samples}

Returns the dictionary of samples of a polynomial system `F`.
"""
samples(F::SampledSystem) = F.samples

"""
    ninstances(F::SampledSystem) -> Int

Returns the number of sampled instances of `F`.
"""
ninstances(F::SampledSystem) = sum([ninstances(s) for s in values(samples(F))])

"""
    nsamples(F::SampledSystem) -> Int

Returns the number of samples of `F`. Notice that `ninstances(F)*nsolutions(F)` doesn't
have to be equal to `nsamples(F)`.
"""
function nsamples(F::SampledSystem)
    return sum([nsamples(s) for s in values(samples(F))])
end

"""
    monodromy_permutations(F::SampledSystem) -> Vector{Vector{Int}}

Returns the vector of monodromy permutations of `F` obtained by [`run_monodromy`](@ref).
"""
monodromy_permutations(F::SampledSystem) = F.mon_info.monodromy_permutations

"""
    block_partitions(F::SampledSystem) -> Vector{Vector{Vector{Int}}}

Returns the vector of all block partitions of the solutions of `F`.
"""
block_partitions(F::SampledSystem) = F.mon_info.block_partitions

"""
    deck_permutations(F::SampledSystem) -> Vector{Vector{Int}}

Returns the vector of deck permutations of the solutions (actions of deck transformations
on the solutions) of `F`.
"""
deck_permutations(F::SampledSystem) = F.mon_info.deck_permutations

(F::SampledSystem)(
    x₀::AbstractVector{<:Number},
    p₀::AbstractVector{<:Number}
) = F.system(x₀, p₀)

function _filter_permutations(perms::Matrix{Int})::Vector{Vector{Int}}
    nsols = length(perms[:,1])
    return filter(
        x->!(0 in x) && (length(unique(x)) == nsols),
        eachcol(perms)
    )
end

function SampledSystem(F::System, MR::MonodromyResult)
    sols, params  = hcat(HC.solutions(MR)...), MR.parameters
    n_sols = size(sols, 2)

    if n_sols == 1
        @warn "Monodromy result has only 1 solution, no monodromy group available"
        return SampledSystem(
            F,
            MonodromyInfo(),
            Dict([1] => Samples(sols, params))
        )
    end

    monodromy_permutations = _filter_permutations(HC.permutations(MR))
    block_partitions = all_block_partitions(to_group(monodromy_permutations))
    deck_permutations = to_permutations(centralizer(monodromy_permutations))

    return SampledSystem(
        F,
        MonodromyInfo(
            n_sols,
            monodromy_permutations,
            block_partitions,
            deck_permutations
        ),
        Dict(Vector(1:n_sols) => Samples(sols, params))
    )
end

function Base.show(io::IO, F::SampledSystem)
    println(io, "SampledSystem with $(phrase(nsamples(F), "sample"))")
    print(io, " $(phrase(nunknowns(F), "unknown")): ", join(unknowns(F), ", "))
    if !isempty(parameters(F))
        print(io, "\n $(phrase(nparameters(F), "parameter")): ", join(parameters(F), ", "))
    end
    print(io, "\n\n")
    println(io, " number of solutions: $(nsolutions(F))")
    print(io, " sampled instances: $(ninstances(F))")
    # print(io, " deck permutations: $(length(deck_permutations(F)))")
end

function random_samples(samples::Samples)
    instance_id = rand(1:ninstances(samples))
    return M2VV(samples.solutions[:, :, instance_id]), samples.parameters[:, instance_id]
end

all_solutions_samples(F::SampledSystem) = samples(F)[Vector(1:nsolutions(F))]

function random_samples(
    F::SampledSystem;
    path_ids::Vector{Int}
)
    samples = all_solutions_samples(F)
    instance_id = rand(1:ninstances(samples))
    return M2VV(samples.solutions[:, path_ids, instance_id]), samples.parameters[:, instance_id]
end


"""
    run_monodromy(F::Union{System, AbstractSystem}, xp₀=nothing; options...) -> SampledSystem

Runs [`monodromy_solve`]($(MONODROMY_SOLVE_REF)) on a given polynomial system `F` with starting
solutions `xp₀[1]` and parameters `xp₀[2]` (if given).

```julia-repl
julia> @var x a b;

julia> F = System([x^3+a*x+b]; variables=[x], parameters=[a,b]);

julia> F = run_monodromy(F, ([[1]], [1,-2]); max_loops_no_progress = 10)
SampledSystem with 3 samples
 1 unknown: x
 2 parameters: a, b

 number of solutions: 3
 sampled instances: 1
```
"""
function run_monodromy(
    F::Union{System, AbstractSystem},
    xp₀::Union{Nothing, Tuple{AbstractVector{<:AbstractVector{<:Number}}, AbstractVector{<:Number}}}=nothing;
    options...
)
    if isnothing(xp₀)
        MR = HC.monodromy_solve(F; permutations=true, options...)
    else
        sols, p₀ = xp₀
        MR = HC.monodromy_solve(F, sols, ComplexF64.(p₀); permutations=true, options...)
    end
    if length(HC.solutions(MR)) == 1
        error("Only 1 solution found, no monodromy group available. Try running again...")
    end
    return SampledSystem(F, MR)
end

"""
    run_monodromy(F::SampledSystem, xp₀=nothing; options...) -> SampledSystem

Reruns [`monodromy_solve`]($(MONODROMY_SOLVE_REF)) on a given sampled polynomial system `F`.
"""
function run_monodromy(
    F::SampledSystem,
    xp₀::Union{Nothing, Tuple{AbstractVector{<:AbstractVector{<:Number}}, AbstractVector{<:Number}}}=nothing;
    options...
)
    if isnothing(xp₀)
        sols, p₀ = random_samples(F; path_ids=Vector(1:nsolutions(F)))
    else
        sols, p₀ = xp₀  # TODO: do we need this?
    end
    MR = HC.monodromy_solve(F.system, sols, ComplexF64.(p₀); permutations=true, options...)
    if length(HC.solutions(MR)) == 1
        error("Only 1 solution found, no monodromy group available. Try running again...")
    end
    return SampledSystem(F.system, MR)
end

function extract_samples(
    results::Vector{Tuple{Result, Vector{ComplexF64}}},
    F::SampledSystem;
    resample::Bool=false
)
    n_tracked = ntracked(results[1][1])
    n_instances = length(results)
    all_sols = zeros(ComplexF64, nunknowns(F), n_tracked, n_instances)
    all_params = zeros(ComplexF64, nparameters(F), n_instances)
    k = 1
    for (res, p) in results
        sols = HC.solutions(res)
        if length(sols) == n_tracked
            all_sols[:, :, k] = hcat(sols...)
            all_params[:, k] = p
            k += 1
        elseif !resample
            error("Number of solutions in the $(k)-th result is $(length(sols)), expected $(n_tracked)")
        end
    end
    for i in k:n_instances
        while true
            instance_id = rand(1:i-1)  # TODO: what if i == 1?
            p₀ = all_params[:, instance_id]
            sols₀ = M2VV(all_sols[:, :, instance_id])
            p₁ = randn(ComplexF64, nparameters(F))
            res = HC.solve(
                F.system,
                sols₀,
                start_parameters = p₀,
                target_parameters = p₁
            )
            sols = HC.solutions(res)
            if length(sols) == n_tracked
                all_sols[:, :, i] = hcat(sols...)
                all_params[:, i] = p₁
                break
            end
        end
    end
    return all_sols, all_params
end

"""
    sample!(F::SampledSystem; path_ids=Vector(1:nsolutions(F)), n_instances=1) -> SampledSystem

Uses [`solve`]($(SOLVE_REF)) method to track the solutions of a poynomial system `F` with ids
defined by `path_ids` to `n_instances` random parameters.
"""
function sample!(
    F::SampledSystem;
    path_ids::AbstractVector{Int}=1:nsolutions(F),
    n_instances::Int=1
)
    (length(path_ids) == 0 || n_instances ≤ 0) && return F
    p₁s = [randn(ComplexF64, nparameters(F)) for _ in 1:n_instances]
    path_ids = sort(Vector{Int}(path_ids))
    samples = get(F.samples, path_ids, nothing)
    if isnothing(samples)
        sols₀, p₀ = random_samples(F; path_ids=path_ids)
        res = HC.solve(
            F.system,
            sols₀,
            start_parameters = p₀,
            target_parameters = p₁s
        )
        sols, params = extract_samples(res, F; resample=true)
        F.samples[path_ids] = Samples(sols, params)
    else
        sols₀, p₀ = random_samples(samples)
        res = HC.solve(
            F.system,
            sols₀,
            start_parameters = p₀,
            target_parameters = p₁s
        )
        sols, params = extract_samples(res, F; resample=true)
        sols = cat(samples.solutions, sols; dims=3)
        params = cat(samples.parameters, params; dims=2)
        F.samples[path_ids] = Samples(sols, params)
    end
    return F
end

function track_parameter_homotopy(
    F::System,
    (x₀, p₀)::NTuple{2, AbstractVector{<:Number}},
    p₁::AbstractVector{<:Number},
    p_inter::AbstractVector{<:Number} # intermediate parameter
)

    H₁ = ParameterHomotopy(F; start_parameters=p₀, target_parameters=p_inter)
    res = track(Tracker(H₁), x₀)
    if !is_success(res)
        @warn "Tracking was not successful: stopped at t = $(res.t)"
    end

    H₂ = ParameterHomotopy(F; start_parameters=p_inter, target_parameters=p₁)
    res = track(Tracker(H₂), solution(res))
    if !is_success(res)
        @warn "Tracking was not successful: stopped at t = $(res.t)"
    end

    return solution(res)
end