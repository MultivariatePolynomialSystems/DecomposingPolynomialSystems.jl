export SampledSystem,
    MonodromyInfo,
    Samples,
    run_monodromy,
    sample_system!,
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

using HomotopyContinuation: Result, MonodromyResult, nsolutions, ntracked
using HomotopyContinuation: ParameterHomotopy, Tracker, track

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
unknowns(F::SampledSystem) = unknowns(F.system)
HC.parameters(F::SampledSystem) = parameters(F.system)
variables(F::System) = vcat(unknowns(F), parameters(F))  # does a different thing than HC.variables
variables(F::SampledSystem) = variables(F.system)

nunknowns(F::SampledSystem) = length(unknowns(F))  # equivalent to HC.nvariables
HC.nparameters(F::SampledSystem) = length(parameters(F))
nvariables(F::SampledSystem) = length(variables(F))  # doesn't extend HC.nvariables, does a different thing

HC.nsolutions(F::SampledSystem) = F.mon_info.n_solutions

function nsamples(F::SampledSystem)
    return sum([nsamples(s) for s in values(samples(F))])
end

function ninstances(F::SampledSystem)
    return sum([ninstances(s) for s in values(samples(F))])
end

samples(F::SampledSystem) = F.samples
monodromy_permutations(F::SampledSystem) = F.mon_info.monodromy_permutations
block_partitions(F::SampledSystem) = F.mon_info.block_partitions
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
    println(io, " sampled instances: $(ninstances(F))")
    print(io, " deck permutations: $(length(deck_permutations(F)))")
end

function random_samples(samples::Samples)
    instance_id = rand(1:ninstances(samples))
    return M2VV(samples.solutions[:, :, instance_id]), samples.parameters[:, instance_id]
end

function random_samples(
    F::SampledSystem;
    path_ids::Vector{Int}
)
    samples = F.samples[Vector(1:nsolutions(F))]
    instance_id = rand(1:ninstances(samples))
    return M2VV(samples.solutions[:, path_ids, instance_id]), samples.parameters[:, instance_id]
end

function run_monodromy(
    F::System,
    xp₀::Union{Nothing, NTuple{2, AbstractVector{<:Number}}}=nothing;
    options...
)
    if isnothing(xp₀)
        MR = HC.monodromy_solve(F; permutations=true, options...)
    else
        x₀, p₀ = xp₀
        MR = HC.monodromy_solve(F, [x₀], p₀; permutations=true, options...)
    end
    if length(HC.solutions(MR)) == 1
        error("Only 1 solution found, no monodromy group available. Try running again...")
    end
    return SampledSystem(F, MR)
end

function run_monodromy(F::SampledSystem; options...)
    xs, p = random_samples(F; path_ids=Vector(1:nsolutions(F)))
    MR = HC.monodromy_solve(F.system, xs, p; permutations=true, options...)
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
            instance_id = rand(1:i-1)
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

function sample_system!(
    F::SampledSystem;
    path_ids::Vector{Int}=Vector(1:nsolutions(F)),
    n_instances::Int=1
)
    p₁s = [randn(ComplexF64, nparameters(F)) for _ in 1:n_instances]
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
    p₁::AbstractVector{<:Number}
)
    H = ParameterHomotopy(F; start_parameters=p₀, target_parameters=p₁) # straight line between p₀ and p₁
    res = track(Tracker(H), x₀)
    if !is_success(res)
        @warn "Tracking was not successful: stopped at t = $(res.t)"
    end
    return res.solution
end