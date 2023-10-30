export SampledSystem,
    VarietySamples,
    run_monodromy,
    sample_system!,
    unknowns,
    parameters,
    variables,
    n_unknowns,
    n_parameters,
    n_variables

using HomotopyContinuation: Result, MonodromyResult, track

# TODO: think about other ways to represent samples
# TODO: make it parametric with T <: Number?
struct VarietySamples
    solutions::Array{CC, 3} # n_unknowns x n_sols x n_instances
    parameters::Array{CC, 2} # n_params x n_instances
end

# TODO
function Base.show(io::IO, samples::VarietySamples)
    println(io, "VarietySamples")
    println(io, " solutions:")
    print(io, " parameters:")
end

mutable struct SampledSystem
    system::System
    samples::VarietySamples
    monodromy_permutations::Vector{Vector{Int}}
    block_partitions::Vector{Vector{Vector{Int}}}
    deck_permutations::Vector{Vector{Int}}
end

function filter_permutations(perms::Matrix{Int})::Vector{Vector{Int}}
    nsols = length(perms[:,1])
    return filter(
        x->!(0 in x) && (length(unique(x)) == nsols),
        eachcol(perms)
    )
end

function SampledSystem(F::System, MR::MonodromyResult)
    sols = hcat(HC.solutions(MR)...)
    sols, params = reshape(sols, size(sols)..., 1), reshape(MR.parameters, :, 1)

    # TODO: throw warning? Error?
    if size(sols, 2) == 1
        return SampledSystem(F,
            VarietySamples(sols, params),
            [], [], []
        )
    end

    monodromy_permutations = filter_permutations(HC.permutations(MR))
    block_partitions = all_block_partitions(to_group(monodromy_permutations))
    deck_permutations = to_permutations(centralizer(monodromy_permutations))

    return SampledSystem(
        F,
        VarietySamples(sols, params),
        monodromy_permutations,
        block_partitions,
        deck_permutations
    )
end

function Base.show(io::IO, F::SampledSystem)
    sols = F.samples.solutions
    n_samples = size(sols, 2)*size(sols, 3)
    println(io, "SampledSystem with $(phrase(n_samples, "sample"))")
    print(io, " $(phrase(n_unknowns(F), "unknown")): ", join(unknowns(F), ", "))
    if !isempty(parameters(F))
        print(io, "\n $(phrase(n_parameters(F), "parameter")): ", join(parameters(F), ", "))
    end
    print(io, "\n\n")
    println(io, " sampled instances: $(size(sols, 3))")
    println(io, " solutions per instance: $(size(sols, 2))")
    print(io, " deck permutations: $(length(F.deck_permutations))")
end

unknowns(F::System) = HC.variables(F)
unknowns(F::SampledSystem) = unknowns(F.system)
HC.parameters(F::SampledSystem) = parameters(F.system)
variables(F::System) = vcat(unknowns(F), parameters(F))
variables(F::SampledSystem) = variables(F.system)

# TODO: remove underscore?
n_unknowns(F::SampledSystem) = length(unknowns(F))
n_parameters(F::SampledSystem) = length(parameters(F))
n_variables(F::SampledSystem) = length(variables(F))

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

# TODO
function run_monodromy(F::SampledSystem; options...)

end

function extract_samples(data_points::Vector{Tuple{Result, Vector{CC}}}, F::SampledSystem)::Tuple{Array{CC, 3}, Array{CC, 2}}
    p0 = F.samples.parameters[:, 1]
    sols0 = M2VV(F.samples.solutions[:, :, 1])
    n_unknowns, n_sols, _ = size(F.samples.solutions)
    n_params = length(p0)
    n_instances = length(data_points)
    all_sols = zeros(CC, n_unknowns, n_sols, n_instances)
    all_params = zeros(CC, n_params, n_instances)
    for i in 1:n_instances
        sols = HC.solutions(data_points[i][1])
        params = data_points[i][2]
        while length(sols) != n_sols
            # println("Tracking solutions again...")
            res = HC.solve(F.system,
                sols0,
                start_parameters = p0,
                target_parameters = [randn(CC, n_params)]
            )
            sols = HC.solutions(res[1][1])
            params = res[1][2]
        end
        all_sols[:, :, i] = hcat(sols...)
        all_params[:, i] = params
    end
    return (all_sols, all_params)
end

function sample_system!(F::SampledSystem, target_params::AbstractVector{<:Number})
    instance_id = rand(1:size(F.samples.solutions, 3))
    p0 = F.samples.parameters[:, instance_id]
    sols = M2VV(F.samples.solutions[:, :, instance_id])
    data_points = HC.solve(
        F.system,
        sols,
        start_parameters = p0,
        target_parameters = [target_params]
    )
    
    n_unknowns, n_sols, _ = size(F.samples.solutions)
    sols = zeros(CC, n_unknowns, n_sols, 1)
    sols[:, :, 1] = hcat(HC.solutions(data_points[1][1])...) # TODO: what if n_sols is different?
    n_params = size(F.samples.parameters, 1)
    params = zeros(CC, n_params, 1)
    params[:, 1] = data_points[1][2]
    
    all_sols = cat(F.samples.solutions, sols, dims=3)
    all_params = cat(F.samples.parameters, params, dims=2)

    F.samples = VarietySamples(all_sols, all_params)
end

# n_instances is the desired number of sampled instances in total
function sample_system!(F::SampledSystem, n_instances::Int)
    n_computed_instances = size(F.samples.parameters, 2)
    if n_computed_instances < n_instances
        p0 = F.samples.parameters[:, 1]
        sols = M2VV(F.samples.solutions[:, :, 1])
        target_params = [randn(CC, length(p0)) for _ in 1:(n_instances-n_computed_instances)]

        # println("Solving ", n_instances-n_computed_instances, " instances by homotopy continutation...")
        data_points = HC.solve(
            F.system,
            sols,
            start_parameters = p0,
            target_parameters = target_params
        )

        solutions, parameters = extract_samples(data_points, F)
        all_sols = cat(F.samples.solutions, solutions, dims=3)
        all_params = cat(F.samples.parameters, parameters, dims=2)
        F.samples = VarietySamples(all_sols, all_params)
    end

    return Vector((n_computed_instances+1):n_instances)  # TODO: remove return?
end

function HC.track(
    (x₀, p₀)::NTuple{2, AbstractVector{<:Number}},
    p₁::AbstractVector{<:Number},
    F::System
)
    data_points = HC.solve(
        F,
        [x₀],
        start_parameters = p₀,
        target_parameters = [p₁]
    )
    return HC.solutions(data_points[1][1])[1]
end