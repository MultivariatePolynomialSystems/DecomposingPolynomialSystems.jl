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

using HomotopyContinuation:
    MonodromyResult,
    Result,
    monodromy_solve,
    solutions

# TODO: think about other ways to represent samples
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

# TODO: Do we need this constructor?
SampledSystem() = SampledSystem(
    System([]),
    VarietySamples(Array{CC}(undef, 0, 0, 0), Array{CC}(undef, 0, 0)),
    Array{Int}(undef, 0, 0),
    [], []
)

function SampledSystem(F::System, MR::MonodromyResult)
    sols, p₀ = HomotopyContinuation.solutions(MR), MR.parameters
    
    solsM = hcat(sols...)
    solutions, parameters = reshape(solsM, size(solsM)..., 1), reshape(p₀, length(p₀), 1)

    # TODO: throw warning? Error?
    if length(sols) == 1
        return SampledSystem(F,
            VarietySamples(solutions, parameters),
            [], [], []
        )
    end

    monodromy_permutations = filter_permutations(HomotopyContinuation.permutations(MR))
    block_partitions = all_block_partitions(permutations_to_group(monodromy_permutations))
    deck_permutations = group_to_permutations(centralizer(monodromy_permutations))

    return SampledSystem(
        F,
        VarietySamples(solutions, parameters),
        monodromy_permutations,
        block_partitions,
        deck_permutations
    )
end

function Base.show(io::IO, F::SampledSystem)
    sols = F.samples.solutions
    n_samples = size(sols, 2)*size(sols, 3)
    println(io, "SampledSystem with $(n_samples) samples")
    print(io, " $(n_unknowns(F)) unknowns: ", join(unknowns(F), ", "))
    if !isempty(parameters(F))
        print(io, "\n $(n_parameters(F)) parameters: ", join(parameters(F), ", "))
    end
    print(io, "\n\n")
    println(io, " number of solutions: $(size(sols, 2))")
    println(io, " number of instances: $(size(sols, 3))")
    print(io, " number of deck permutations: $(length(F.deck_permutations))")
end

unknowns(F::SampledSystem) = F.system.variables
parameters(F::SampledSystem) = F.system.parameters
variables(F::SampledSystem) = vcat(unknowns(F), parameters(F))

# TODO: remove underscore?
n_unknowns(F::SampledSystem) = length(unknowns(F))
n_parameters(F::SampledSystem) = length(parameters(F))
n_variables(F::SampledSystem) = length(variables(F))

function run_monodromy(F::System; options...)::SampledSystem
    MR = monodromy_solve(F; permutations=true, options...)
    if length(solutions(MR)) == 1
        error("Just one solution was found, no monodromy group available. Try running again...")
    end
    return SampledSystem(F, MR)
end

function run_monodromy(F::System, (x₀, p₀)::Tuple{Vector{CC}, Vector{CC}}; options...)::SampledSystem
    MR = monodromy_solve(F, [x₀], p₀; permutations=true, options...)
    if length(solutions(MR)) == 1
        error("No additional solutions were found, no monodromy group available. Try running again...")
    end
    return SampledSystem(F, MR)
end

# TODO
function run_monodromy(F::SampledSystem; options...)::SampledSystem

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
        sols = HomotopyContinuation.solutions(data_points[i][1])
        params = data_points[i][2]
        while length(sols) != n_sols
            # println("Tracking solutions again...")
            res = HomotopyContinuation.solve(F.system,
                sols0,
                start_parameters = p0,
                target_parameters = [randn(CC, n_params)]
            )
            sols = HomotopyContinuation.solutions(res[1][1])
            params = res[1][2]
        end
        all_sols[:, :, i] = hcat(sols...)
        all_params[:, i] = params
    end
    return (all_sols, all_params)
end

function sample_system!(F::SampledSystem, target_params::Vector{CC})
    instance_id = rand(1:size(F.samples.solutions, 3))
    p0 = F.samples.parameters[:, instance_id]
    sols = M2VV(F.samples.solutions[:, :, instance_id])
    data_points = HomotopyContinuation.solve(
        F.system,
        sols,
        start_parameters = p0,
        target_parameters = [target_params]
    )
    
    n_unknowns, n_sols, _ = size(F.samples.solutions)
    solutions = zeros(CC, n_unknowns, n_sols, 1)
    solutions[:, :, 1] = hcat(HomotopyContinuation.solutions(data_points[1][1])...) # TODO: what if n_sols is different?
    n_params = size(F.samples.parameters, 1)
    parameters = zeros(CC, n_params, 1)
    parameters[:, 1] = data_points[1][2]
    
    all_sols = cat(F.samples.solutions, solutions, dims=3)
    all_params = cat(F.samples.parameters, parameters, dims=2)

    F.samples = VarietySamples(all_sols, all_params)
end

# n_instances is the desired number of sampled instances in total
function sample_system!(F::SampledSystem, n_instances::Int)
    n_computed_instances = size(F.samples.parameters, 2)
    if n_computed_instances < n_instances
        p0 = F.samples.parameters[:, 1]
        sols = M2VV(F.samples.solutions[:, :, 1])
        target_params = [randn(CC, length(p0)) for _ in 1:(n_instances-n_computed_instances)]

        println("Solving ", n_instances-n_computed_instances, " instances by homotopy continutation...")
        data_points = HomotopyContinuation.solve(
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