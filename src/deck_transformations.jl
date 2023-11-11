export Tolerances,
    DeckTransformation,
    DeckTransformationGroup,
    symmetries_fixing_parameters_dense!,
    symmetries_fixing_parameters_graded!,
    symmetries_fixing_parameters!,
    symmetries_fixing_parameters,
    _deck_action, _deck_commutes_with_scaling,
    _scalings_commuting_with_deck

@kwdef struct Tolerances
    nullspace_atol::Float64=0
    nullspace_rtol::Float64=0
    rref_tol::Float64=1e-5
    sparsify_tol::Float64=1e-5
end

MiExpression = Union{Missing, Expression}

struct DeckTransformation
    exprs::Vector{MiExpression}
    unknowns::Vector{Variable}
    parameters::Vector{Variable}

    function DeckTransformation(exprs, unknowns, parameters)
        # TODO: verify args
        return new(exprs, unknowns, parameters)
    end
end

Base.getindex(dt::DeckTransformation, inds...) = getindex(dt.exprs, inds...)

function Base.show(io::IO, dt::DeckTransformation)
    println(
        io,
        "DeckTransformation: acts on $(phrase(length(dt.unknowns), "unknown")),",
        " fixes $(phrase(length(dt.parameters), "parameter"))",
    )
    println(io, " action:")
    for i in 1:length(dt.exprs)
        print(io, "  ", dt.unknowns[i], " ↦ ", dt.exprs[i])
        i < length(dt.exprs) && print(io, "\n")
    end
end

struct DeckTransformationGroup
    maps::Vector{DeckTransformation}
    group::GapObj
    F::SampledSystem
end

function DeckTransformationGroup(F::SampledSystem)
    symmetries = _init_symmetries(length(F.deck_permutations), unknowns(F))
    return DeckTransformationGroup(symmetries, F)
end

function DeckTransformationGroup(
    symmetries::Vector{Vector{MiExpression}},
    F::SampledSystem
)
    action = [DeckTransformation(symmetry, unknowns(F), parameters(F)) for symmetry in symmetries]
    return DeckTransformationGroup(action, to_group(F.deck_permutations), F)
end

function Base.show(io::IO, deck::DeckTransformationGroup)
    println(io, "DeckTransformationGroup of order $(length(deck.maps))")
    println(io, " structure: ", order(deck.group) == 1 ? "trivial" : group_structure(deck.group))
    print(io, " action:")
    for i in eachindex(deck.maps)
        println(io, "\n  ", to_ordinal(i), " map:")
        for (j, var) in enumerate(unknowns(deck.F))  # action on parameters is trivial, don't show it
            print(io, "   ", var, " ↦ ", deck.maps[i][j])
            j < length(unknowns(deck.F)) && print(io, "\n")
        end
    end
end

Base.getindex(deck::DeckTransformationGroup, inds...) = getindex(deck.maps, inds...)

# TODO: extend to scalings::ScalingGroup?
function _denom_deg(num_deg::Vector{Int}, grading::Grading, var_id::Int)
    denom_deg = zeros(Int, length(num_deg))
    U₀ = grading.free_part
    if !isnothing(U₀)
        k = size(U₀, 1)
        denom_deg[1:k] = num_deg[1:k] - U₀[:, var_id]
    end
    for (sᵢ, Uᵢ) in grading.mod_part
        n_scalings = size(Uᵢ, 1)
        denom_deg[k+1:k+n_scalings] = mod.(num_deg[k+1:k+n_scalings] - Uᵢ[:, var_id], sᵢ)
        k += n_scalings
    end
    return denom_deg
end

function _deck_vandermonde_matrix(
    deck_permutation::Vector{Int},
    function_id::Int,
    solutions::AbstractArray{T, 3},
    eval_num_mons::AbstractArray{T, 3},
    eval_denom_mons::AbstractArray{T, 3}
) where {T<:Complex}

    _, n_sols, n_instances = size(solutions)
    n_num_mons = size(eval_num_mons, 1)
    n_denom_mons = size(eval_denom_mons, 1)

    A = zeros(T, n_instances*n_sols, n_num_mons+n_denom_mons)
    @assert size(A, 1) >= size(A, 2)

    for i in 1:n_sols
        v = solutions[function_id, deck_permutation[i], :]
        rows = ((i-1)*n_instances+1):(i*n_instances)
        A[rows, 1:n_num_mons] = transpose(eval_num_mons[:, i, :])
        A[rows, (n_num_mons+1):end] = -transpose(eval_denom_mons[:, i, :]).*v
    end
    return A
end

function _all_interpolated(symmetries::Vector{Vector{MiExpression}})
    all_interpolated = true
    for symmetry in symmetries
        for expr in symmetry
            if ismissing(expr)
                all_interpolated = false
                break
            end
        end
    end
    return all_interpolated
end

function _init_symmetries(n_symmetries::Int, unknowns::Vector{Variable})
    symmetries = [[missing for j in eachindex(unknowns)] for i in 1:n_symmetries]
    symmetries = Vector{Vector{MiExpression}}(symmetries)
    symmetries[1] = Expression.(unknowns)  # set the first to the identity
    return symmetries
end

function _interpolate_deck_function(
    deck_permutation::Vector{Int},
    function_id::Int,
    solutions::AbstractArray{T, 3},
    eval_num_mons::AbstractArray{T, 3},
    eval_denom_mons::AbstractArray{T, 3},
    num_mons::MonomialVector,
    denom_mons::MonomialVector,
    tols::Tolerances;
    logging::Bool=false
) where {T<:Complex}

    logging && println(
        "Creating vandermonde matrix of size ",
        (prod([size(solutions,2), size(solutions,3)]), length(num_mons)+length(denom_mons))
    )
    A = _deck_vandermonde_matrix(
        deck_permutation,
        function_id,
        solutions,
        eval_num_mons,
        eval_denom_mons
    )

    logging && println("Computing nullspace...")
    if tols.nullspace_rtol == 0
        N = nullspace(A, atol=tols.nullspace_atol)
    else
        N = nullspace(A, atol=tols.nullspace_atol, rtol=tols.nullspace_rtol)
    end
    coeffs = transpose(N)
    logging && println("Size of the transposed nullspace: ", size(coeffs))

    if size(coeffs, 1) == 0 return missing end

    logging && println("Computing the reduced row echelon form of the transposed nullspace...\n")
    coeffs = rref(coeffs, tols.rref_tol)
    
    sparsify!(coeffs, tols.sparsify_tol; digits=1)
    coeffs = remove_zero_nums_and_denoms(coeffs, num_mons, denom_mons)
    if size(coeffs, 1) == 0 return missing end

    coeffs = good_representative(coeffs)
    return rational_function(coeffs, num_mons, denom_mons; tol=tols.sparsify_tol)
end

function _interpolate_deck_function(
    deck_permutation::Vector{Int},
    function_id::Int,
    solutions::AbstractArray{T, 3},
    eval_mons::AbstractArray{T, 3},
    mons::MonomialVector,
    tols::Tolerances;
    logging::Bool=false
) where {T<:Complex}

    return _interpolate_deck_function(
        deck_permutation,
        function_id,
        solutions,
        eval_mons,
        eval_mons,
        mons,
        mons,
        tols;
        logging=logging
    )
end

function symmetries_fixing_parameters_graded!(
    F::SampledSystem,
    scalings::ScalingGroup,
    mons::MonomialVector,
    classes::Dict{Vector{Int}, Vector{Int}};
    tols::Tolerances=Tolerances(),
    logging::Bool=false
)
    
    max_n_mons = max(length.(collect(values(classes)))...)  # size of the largest class
    n_unknowns, n_sols, _ = size(F.samples.solutions)  # TODO: what if n_sols is huge?
    n_instances = Int(ceil(2/n_sols*max_n_mons))

    C = F.deck_permutations
    symmetries = _init_symmetries(length(C), unknowns(F))

    sample_system!(F, n_instances)
    
    for (num_deg, num_ids) in classes
        num_mons = mons[num_ids]
        eval_num_mons = nothing
        for i in 1:n_unknowns
            denom_deg = _denom_deg(num_deg, scalings.grading, i)  # i-th variable
            denom_ids = get(classes, denom_deg, nothing)
            if !isnothing(denom_ids)
                denom_mons = mons[denom_ids]
                g = gcd(vcat(num_mons, denom_mons))
                if isone(g) && !only_param_dep(vcat(num_mons, denom_mons), Vector(1:n_unknowns))
                    if isnothing(eval_num_mons)
                        eval_num_mons = HC.evaluate(num_mons, F.samples)
                    end
                    eval_denom_mons = HC.evaluate(denom_mons, F.samples)
                    for (j, symmetry) in enumerate(symmetries)
                        if ismissing(symmetry[i])
                            symmetry[i] = _interpolate_deck_function(
                                C[j],
                                i,
                                F.samples.solutions,
                                eval_num_mons,
                                eval_denom_mons,
                                num_mons,
                                denom_mons,
                                tols;
                                logging=logging
                            )
                            if logging && !ismissing(symmetry[i])
                                printstyled(
                                    "Good representative for the ",
                                    to_ordinal(j),
                                    " symmetry, variable ",
                                    unknowns(F)[i],
                                    ":\n",
                                    color=:red
                                )
                                println(symmetry[i])
                            end
                        end
                    end
                end
            end
        end

        if _all_interpolated(symmetries)
            logging && printstyled("--- All symmetries are interpolated ---\n", color=:blue)
            return DeckTransformationGroup(symmetries, F)
        end
    end

    return DeckTransformationGroup(symmetries, F)
end

function symmetries_fixing_parameters_graded!(
    F::SampledSystem,
    scalings::ScalingGroup;
    degree_bound::Integer=1,
    tols::Tolerances=Tolerances(),
    logging::Bool=false
)

    mons = MonomialVector{Int8}(scalings.vars; degree=degree_bound)
    classes = to_classes(mons, scalings.grading)
    return symmetries_fixing_parameters_graded!(
        F,
        scalings,
        mons,
        classes;
        tols=tols,
        logging=logging
    )
end

function symmetries_fixing_parameters_dense!(
    F::SampledSystem; 
    degree_bound::Integer=1,
    param_dep::Bool=true,
    tols::Tolerances=Tolerances(),
    logging::Bool=false
)

    n_unknowns, n_sols, _ = size(F.samples.solutions)  # TODO: what if n_sols is huge?
    vars = param_dep ? variables(F) : unknowns(F)  # TODO: rename vars --> interp_vars?

    C = F.deck_permutations
    symmetries = _init_symmetries(length(C), unknowns(F))

    for d in 1:degree_bound
        logging && printstyled("Started interpolation for degree = ", d, "...\n"; color=:green)
        mons = MonomialVector{Int8}(vars, degree=d)
        n_instances = Int(ceil(2/n_sols*length(mons)))
        sample_system!(F, n_instances)

        logging && println("Evaluating monomials...\n")
        evaluated_mons = HC.evaluate(mons, F.samples)

        for (i, symmetry) in enumerate(symmetries)
            logging && printstyled("Interpolating the ", i, "-th symmetry map...\n"; color=:blue)
            for j in 1:n_unknowns
                if ismissing(symmetry[j])
                    symmetry[j] = _interpolate_deck_function(
                        C[i],
                        j,
                        F.samples.solutions,
                        evaluated_mons,
                        mons,
                        tols;
                        logging=logging
                    )
                    if logging && !ismissing(symmetry[j])
                        printstyled(
                            "Good representative for the ",
                            i,
                            "-th symmetry, variable ",
                            unknowns(F)[j],
                            ":\n";
                            color=:red
                        )
                        println(symmetry[j])
                    end
                end
            end
        end
    
        if _all_interpolated(symmetries)
            logging && printstyled("--- All symmetries are interpolated ---\n"; color=:blue)
            return DeckTransformationGroup(symmetries, F)
        end
    end

    return DeckTransformationGroup(symmetries, F)
end

to_CC(scaling::Tuple{Int, Vector{Int}}) = [cis(2*pi*k/scaling[1]) for k in scaling[2]]

function _deck_action(
    deck_permutation::Vector{Int},
    (x₀, p₀)::NTuple{2, AbstractVector{<:Number}},
    F::SampledSystem;
    tol::Real=1e-5
)
    sols, params = F.samples.solutions, F.samples.parameters
    instance_id = rand(1:size(sols, 3))
    p₁ = params[:, instance_id]
    x₁ = track_parameter_homotopy(F.system, (x₀, p₀), p₁)  # along path γ in the parameter space
    x₁_id = findfirst(x₁, sols[:,:,instance_id]; tol=tol)
    isnothing(x₁_id) && return nothing
    Ψ_x₁ = sols[:,deck_permutation[x₁_id], instance_id]
    Ψ_x₀ = track_parameter_homotopy(F.system, (Ψ_x₁, p₁), p₀)  # should be along the same path γ
    return Ψ_x₀
end

# supposes scaling is a symmetry of F
function _deck_commutes_with_scaling(
    deck_permutation::Vector{Int},
    scaling::Tuple{Int, Vector{Int}},
    F::SampledSystem;
    tol::Real=1e-5
)
    sols, params = F.samples.solutions, F.samples.parameters
    inst_id = rand(1:size(sols, 3))
    p₀ = params[:, inst_id]
    Φ_p₀ = to_CC(scaling)[end-n_parameters(F)+1:end].*p₀
    sol_id = rand(1:size(sols, 2))
    x₀ = sols[:, sol_id, inst_id]
    Ψ_x₀ = sols[:, deck_permutation[sol_id], inst_id]
    Φ_x₀ = to_CC(scaling)[1:n_unknowns(F)].*x₀
    ΦΨ_x₀ = to_CC(scaling)[1:n_unknowns(F)].*Ψ_x₀
    ΨΦ_x₀ = _deck_action(deck_permutation, (Φ_x₀, Φ_p₀), F; tol=tol)
    isnothing(ΨΦ_x₀) && return false
    return norm(ΦΨ_x₀-ΨΦ_x₀)<tol
end

function _all_deck_commute(
    F::SampledSystem,
    scaling::Tuple{Int, Vector{Int}};
    tol::Real=1e-5
)
    all_commute = true
    for deck_permutation in F.deck_permutations
        if !_deck_commutes_with_scaling(deck_permutation, scaling, F; tol=tol)
            all_commute = false
            break
        end
    end
    return all_commute
end

function _scalings_commuting_with_deck(F::SampledSystem, scalings::ScalingGroup)
    grading = scalings.grading
    final_grading = Grading(grading.free_part, [])
    for (sᵢ, Uᵢ) in grading.mod_part
        Vᵢ = Array{Int}(undef, 0, size(Uᵢ, 2))
        # TODO: Uᵢ ↦ all linear combinations of rows of Uᵢ (might not commute with 2 gens, but commutes with their combination)
        for j in axes(Uᵢ, 1)
            if _all_deck_commute(F, (sᵢ, Uᵢ[j, :]))
                Vᵢ = [Vᵢ; hcat(Uᵢ[j, :]...)]
            end
        end
        if size(Vᵢ, 1) > 0
            push!(final_grading.mod_part, (sᵢ, Vᵢ))
        end
    end
    return ScalingGroup(reduce(final_grading), scalings.vars)
end

function symmetries_fixing_parameters!(
    F::SampledSystem;
    degree_bound::Integer=1,
    param_dep::Bool=true,
    tols::Tolerances=Tolerances(),
    logging::Bool=false
)

    if length(F.deck_permutations) == 1 # trivial group of symmetries
        return DeckTransformationGroup(F)
    end

    scalings = _scalings_commuting_with_deck(F, scaling_symmetries(F))
    scalings = param_dep ? scalings : restrict_scalings(scalings, unknowns(F))  # TODO: justify!
    if isempty(scalings.grading)
        logging && printstyled("Running dense version...\n", color=:green)
        return symmetries_fixing_parameters_dense!(
            F;
            degree_bound=degree_bound,
            param_dep=param_dep,
            tols=tols,
            logging=logging
        )
    else
        logging && printstyled("Running graded version...\n", color=:green)
        return symmetries_fixing_parameters_graded!(
            F,
            scalings;
            degree_bound=degree_bound,
            tols=tols,
            logging=logging
        )
    end
end

"""
    symmetries_fixing_parameters(F::System; degree_bound=1, param_dep=true)

Given a polynomial system F returns the group of symmetries 
of `F` that fix the parameters. The keyword
argument `degree_bound` is used to set the upper bound for the
degrees of numerator and denominator polynomials in expressions
for the symmetries. The `param_dep` keyword argument specifies
whether to consider functions of the symmetries to be dependent
on the parameters of `F`.

```julia-repl
julia> @var x[1:2] p[1:2];

julia> F = System([x[1]^2 - x[2]^2 - p[1], 2*x[1]*x[2] - p[2]]; variables=x, parameters=p);

julia> symmetries_fixing_parameters(F; degree_bound=1, param_dep=false)
DeckTransformationGroup of order 4
 structure: C2 x C2
 action:
  1st map:
   x₁ ↦ x₁
   x₂ ↦ x₂
  2nd map:
   x₁ ↦ (0.0 + 1.0*im)*x₂
   x₂ ↦ (0.0 - 1.0*im)*x₁
  3rd map:
   x₁ ↦ (-1.0 + 0.0*im)*x₁
   x₂ ↦ (-1.0 + 0.0*im)*x₂
  4th map:
   x₁ ↦ (0.0 - 1.0*im)*x₂
   x₂ ↦ (0.0 + 1.0*im)*x₁
```
"""
function symmetries_fixing_parameters(  # TODO: extend to take an expression map
    F::System;
    xp₀::Union{Nothing, NTuple{2, AbstractVector{<:Number}}}=nothing,
    degree_bound::Integer=1,
    param_dep::Bool=true,
    tols::Tolerances=Tolerances(),
    monodromy_options::Tuple=(),
    logging::Bool=false
)

    F = run_monodromy(F, xp₀; monodromy_options...)
    return symmetries_fixing_parameters!(
        F;
        degree_bound=degree_bound,
        param_dep=param_dep,
        tols=tols,
        logging=logging
    )
end

