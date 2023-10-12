export RationalMap, SampledSystem, FactorizingMap
export run_monodromy, sample_system!
export evaluate_monomials_at_samples, evaluate_monomials_at_samples_
export unknowns, parameters, variables, n_unknowns, n_parameters, n_variables

Multidegree = Vector{Int8} # TODO: is it OK to suppose degrees are < 2^7 = 128?
Grading = Vector{Tuple{Int, Matrix{Int}}}
NoExpression = Union{Nothing, Expression}

struct RationalMap
    domain_vars::Vector{Variable}
    image_vars::Vector{Variable}
    funcs::Vector{NoExpression}

    function RationalMap(domain_vars, image_vars, funcs)
        # TODO: exclude empty vectors, repetitions in vars
        return new(domain_vars, image_vars, funcs)
    end
end

# TODO: what if vars not in funcs? What if funcs has variables not present in vars?
function RationalMap(vars::Vector{Variable}, funcs::Vector{NoExpression})
    @assert length(vars) == length(funcs) "#vars ≂̸ #funcs, specify image variables"
    return RationalMap(vars, vars, funcs)
end

function Base.show(io::IO, map::RationalMap)
    println(io, "RationalMap: ℂ$(superscriptnumber(length(map.domain_vars))) ⊃ X - - > ℂ$(superscriptnumber(length(map.funcs)))")
    println(io, " action:")
    if map.domain_vars == map.image_vars
        for (i, var) in enumerate(map.domain_vars)
            print(io, "  ", var, " ↦ ", map.funcs[i])
            i < length(map.domain_vars) && print(io, "\n")
        end
    else
        # TODO
    end
end

Base.getindex(map::RationalMap, i::Int) = (map.image_vars[i], map.funcs[i])
function Base.getindex(map::RationalMap, var::Variable)
    id = findfirst(x->x==var, map.image_vars)
    if isnothing(id)
        error("The variable $(var) isn't present in the image variables")
    end
    return map.funcs[id]
end

struct MonomialVector
    mds::Vector{Multidegree}
    vars::Vector{Variable}
end

Base.length(mons::MonomialVector) = length(mons.mds)

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
    
    solsM = VV2M(sols)
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


mutable struct FactorizingMap
    map::Vector{Expression}
    domain::Ref{SampledSystem}
    image::Ref{SampledSystem}
    monodromy_group::GapObj
    # deck_transformations::Vector{Vector{Expression}}
end

function FactorizingMap(map::Vector{Expression})
    return FactorizingMap(map,
        Ref{SampledSystem}(),
        Ref{SampledSystem}(),
        GAP.evalstr( "Group(())" )
    )
end

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
        all_sols[:, :, i] = VV2M(sols)
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
    solutions[:, :, 1] = VV2M(HomotopyContinuation.solutions(data_points[1][1])) # TODO: what if n_sols is different?
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

    return Vector((n_computed_instances+1):n_instances)
end

# supposes each md in mds is a multidegree in both unknowns and parameters
# TODO: The implementation below (with _) is more efficient (approx 2x),
# TODO: since it exploits the sparsity of multidegrees. REMOVE THIS METHOD?
function evaluate_monomials_at_samples(mons::MonomialVector, samples::VarietySamples)::Array{CC, 3}
    solutions = samples.solutions
    parameters = samples.parameters

    n_unknowns, n_sols, n_instances = size(solutions)
    mds = mons.mds
    n_mds = length(mds)

    evaluated_mons = zeros(CC, n_mds, n_sols, n_instances)
    for i in 1:n_instances
        params = parameters[:, i]
        params_eval = [prod(params.^md[n_unknowns+1:end]) for md in mds]
        sols = solutions[:, :, i]
        for j in 1:n_mds
            evaluated_mons[j, :, i] = M2V(prod(sols.^mds[j][1:n_unknowns], dims=1)).*params_eval[j]
        end
    end
    return evaluated_mons
end

# TODO: consider view for slices
function evaluate_monomials_at_samples_(mons::MonomialVector, samples::VarietySamples)::Array{CC, 3}
    solutions = samples.solutions
    parameters = samples.parameters

    n_unknowns, n_sols, n_instances = size(solutions)
    mds = mons.mds
    n_mds = length(mds)

    nonzero_ids_unknowns = [findall(!iszero, md[1:n_unknowns]) for md in mds]
    nonzero_ids_params = [findall(!iszero, md[n_unknowns+1:end]) for md in mds]

    evaluated_mons = zeros(CC, n_mds, n_sols, n_instances)
    for i in 1:n_instances
        params = parameters[:, i]
        sols = solutions[:, :, i]
        for (j, md) in enumerate(mds)
            params_part = params[nonzero_ids_params[j]]
            md_params_part = md[n_unknowns+1:end][nonzero_ids_params[j]]
            params_eval = prod(params_part.^md_params_part)
            sols_part = sols[nonzero_ids_unknowns[j],:]
            md_sols_part = md[1:n_unknowns][nonzero_ids_unknowns[j]]
            evaluated_mons[j, :, i] = M2V(prod(sols_part.^md_sols_part, dims=1)).*params_eval
        end
    end
    return evaluated_mons
end

function multiplication_matrix(F::SampledSystem, f::Expression, B::Vector{Expression}, instance_id::Int)
    sols = F.samples.solutions[:, :, instance_id]
    n_sols = size(sols, 2)
    vars = variables(F.system)
    A = zeros(CC, n_sols, n_sols)
    for i in 1:n_sols
        A[:, i] = express_in_basis(F, f*B[i], B, instance_id=instance_id)
    end
    return sparsify(A, 1e-5)
end

function multiplication_matrix(F::SampledSystem, f::Expression, B::Vector{Expression}; degree::Int)::Matrix{Union{Nothing, Expression}}
    n_sols = size(F.samples.solutions, 2)
    M = Matrix{Union{Nothing, Expression}}(VV2M([[nothing for i in 1:n_sols] for j in 1:n_sols]))
    for i in 1:n_sols
        M[:, i] = express_in_basis(F, f*B[i], B, degree)
    end
    return M
end

function monomial_basis(F::SampledSystem)::Vector{Expression}
    n_sols = size(F.samples.solutions, 2)
    sols = F.samples.solutions[:,:,1]  # solutions for the 1st instance
    unknowns = variables(F.system)
    A = zeros(CC, n_sols, n_sols)
    n_indep = 1
    indep_mons = Vector{Expression}([])

    for i in 0:n_sols-1
        mons = get_monomials_fixed_degree(unknowns, i)
        for j in eachindex(mons)
            A[:, n_indep] = [subs(mons[j], unknowns=>sols[:,k]) for k in 1:n_sols]
            if rank(A, atol=1e-8) == n_indep
                push!(indep_mons, mons[j])
                n_indep += 1
                if n_indep > n_sols
                    return indep_mons
                end
            end
        end
    end
    @warn "Returned monomials don't form a basis!"
    return indep_mons
end

function eval_at_sols(F::SampledSystem, G::Vector{Expression})::Matrix{CC}
    sols = F.samples.solutions[:, :, 1]
    params = F.samples.parameters[:, 1]
    n_sols = size(sols, 2)
    n_elems = length(G)
    vars = vcat(variables(F.system), parameters(F.system))
    A = zeros(CC, n_sols, n_elems)
    for i in 1:n_sols
        A[i, :] = subs(G, vars => vcat(sols[:,i], params))
    end
    return A
end

function is_basis(F::SampledSystem, B::Vector{Expression})::Bool
    sols = F.samples.solutions[:, :, 1]
    n_sols = size(sols, 2)
    if length(B) < n_sols || length(B) > n_sols
        return false
    end
    A = eval_at_sols(F, B)
    N = nullspace(A)
    println("dim null = ", size(N, 2))
    return size(N, 2) == 0
end

function are_LI(F::SampledSystem, G::Vector{Expression})::Bool
    A = eval_at_sols(F, G)
    N = nullspace(A)
    println("dim null = ", size(N, 2))
    return size(N, 2) == 0
end

function express_in_basis(F::SampledSystem, f::Expression, B::Vector{Expression}; instance_id::Int)::Vector{CC}
    sols = F.samples.solutions[:, :, instance_id]
    params = F.samples.parameters[:, instance_id]
    n_sols = size(sols, 2)
    vars = vcat(variables(F.system), parameters(F.system))
    A = zeros(CC, n_sols, n_sols+1)
    for i in 1:n_sols
        A[i, 1:n_sols] = subs(B, vars => vcat(sols[:, i], params))
        A[i, n_sols+1] = -subs(f, vars => vcat(sols[:, i], params))
    end
    c = nullspace(A)
    @assert size(c, 2) == 1
    @assert abs(c[n_sols+1, 1]) > 1e-10
    return sparsify(M2V(p2a(c)), 1e-5)
end

function express_in_basis(F::SampledSystem, f::Expression, B::Vector{Expression}, degree::Int)::Vector{Union{Nothing,Expression}}
    _, n_sols, n_instances = size(F.samples.solutions)
    evals = zeros(CC, n_sols, n_instances)
    for i in 1:n_instances
        evals[:, i] = express_in_basis(F, f, B; instance_id=i)
    end

    params = parameters(F.system)
    mons = get_monomials(params, degree)
    println("n_mons = ", length(mons))
    evaluated_mons = evaluate_monomials_at_samples(mons, F.samples.parameters, params)

    coeffs = Vector{Union{Nothing, Expression}}([nothing for i in 1:n_sols])
    for i in 1:n_sols
        coeffs[i] = interpolate_dense(view(evals, i, :), mons, evaluated_mons, func_type="rational", tol=1e-5)
    end
    return coeffs
end