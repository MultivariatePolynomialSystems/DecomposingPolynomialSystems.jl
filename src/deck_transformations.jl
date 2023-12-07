export Tolerances,
    DeckTransformation,
    DeckTransformationGroup,
    symmetries_fixing_parameters_dense!,
    symmetries_fixing_parameters_graded!,
    symmetries_fixing_parameters!,
    symmetries_fixing_parameters,
    _deck_action, _deck_commutes_with_scaling,
    _scalings_commuting_with_deck,
    _sample_for_deck_computation!,
    _deck_vandermonde_matrix

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

"""
    DeckTransformationGroup

A `DeckTransformationGroup` is the result of deck transformations computation.
"""
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
    return DeckTransformationGroup(action, to_group(deck_permutations(F)), F)
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

function _denom_deg(
    num_deg::SparseVector{Tv,Ti},
    grading::Grading{Tv,Ti},
    var_id::Int
) where {Tv<:Integer,Ti<:Integer}
    denom_deg = spzeros(Tv, Ti, length(num_deg))
    U₀ = grading.free_part
    if !isnothing(U₀)
        denom_deg[1:size(U₀,1)] = num_deg[1:size(U₀,1)] - U₀[:, var_id]
    end
    k = nfree(grading)
    for (sᵢ, Uᵢ) in grading.mod_part
        n_scalings = size(Uᵢ,1)
        denom_deg[(k+1):(k+n_scalings)] = mod.(num_deg[(k+1):(k+n_scalings)] - Uᵢ[:, var_id], sᵢ)
        k += n_scalings
    end
    return denom_deg
end

function to_nconstraints_ninstances(
    path_ids::Dict{Vector{Int}, Vector{Int}},
    samples::Dict{Vector{Int}, Samples}
)
    c_i_dict = Dict{Vector{Int}, NTuple{2, Int}}()
    for key in keys(path_ids)
        c_i_dict[key] = (length(path_ids[key]), ninstances(samples[key]))
    end
    return c_i_dict
end

function total_nconst_ninst(
    c_i_dict::Dict{Vector{Int}, NTuple{2, Int}}
)
    cₜ, iₜ = 0, 0
    for (c, i) in values(c_i_dict)
        cₜ += c*i
        iₜ += i
    end
    return cₜ, iₜ
end

function pick_nconstraints_ninstances(
    c_i_dict::Dict{Vector{Int}, NTuple{2, Int}},
    cₘᵢₙ::Int,
    iₘᵢₙ::Int
)
    picked_c_i = Dict{Vector{Int}, NTuple{2, Int}}()
    cₜ, iₜ = total_nconst_ninst(c_i_dict)
    for (key, (cₖ, iₖ)) in c_i_dict
        cₜ, iₜ = cₜ - cₖ*iₖ, iₜ - iₖ
        Δc, Δi = cₘᵢₙ - cₜ, iₘᵢₙ - iₜ
        if Δc ≤ 0 && Δi ≤ 0
            picked_c_i[key] = (0, 0)
            continue
        end
        i = max(Δi, div(Δc, cₖ, RoundUp))
        i > iₖ && error("Not enough instances or constraints")
        c = max(div(Δc, i, RoundUp), 1)
        c > cₖ && error("Not enough instantces or constraints")
        picked_c_i[key] = (c, i)
        cₘᵢₙ, iₘᵢₙ = cₘᵢₙ - c*i, iₘᵢₙ - i
    end
    return picked_c_i
end

function _deck_vandermonde_matrix(
    deck_permutation::Vector{Int},
    function_id::Int,
    samples::Dict{Vector{Int}, Samples},
    path_ids_for_deck::Dict{Vector{Int}, Vector{Int}},
    eval_num_mons::Dict{Samples, EvaluatedMonomials},
    eval_denom_mons::Dict{Samples, EvaluatedMonomials},
    min_nconstraints::Int,
    min_ninstances::Int
)
    n_num_mons = nmonomials(first(values(eval_num_mons)))
    n_denom_mons = nmonomials(first(values(eval_denom_mons)))
    n_mons = n_num_mons + n_denom_mons
    c_i_dict = to_nconstraints_ninstances(path_ids_for_deck, samples)
    picked_c_i = pick_nconstraints_ninstances(c_i_dict, min_nconstraints, min_ninstances)
    ntotal_constraints, _ = total_nconst_ninst(picked_c_i)
    A = zeros(ComplexF64, ntotal_constraints, n_mons)
    @assert size(A, 1) >= size(A, 2)

    row_offset = 0
    for (sampled_path_ids, samplesₖ) in samples
        n_constraints, n_instances = picked_c_i[sampled_path_ids]
        n_constraints == 0 && continue

        sols = samplesₖ.solutions
        ids_rel = Dict(zip(sampled_path_ids, 1:length(sampled_path_ids)))
        deck_path_ids = path_ids_for_deck[sampled_path_ids]
        eval_num = eval_num_mons[samplesₖ]
        eval_denom = eval_denom_mons[samplesₖ]
        
        for i in 1:n_constraints
            path_id = deck_path_ids[i]  # length(deck_path_ids) ≤ n_constraints always
            deck_values = sols[function_id, ids_rel[deck_permutation[path_id]], 1:n_instances]
            rows = (row_offset+1):(row_offset+n_instances)
            A[rows, 1:nunkn_dep(eval_num)] = transpose(eval_num.unkn_dep[:, ids_rel[path_id], 1:n_instances])
            A[rows, (nunkn_dep(eval_num)+1):n_num_mons] = transpose(eval_num.param_only[:, 1:n_instances])
            A[rows, (n_num_mons+1):(n_num_mons+nunkn_dep(eval_denom))] = -transpose(eval_denom.unkn_dep[:, ids_rel[path_id], 1:n_instances]).*deck_values
            A[rows, (n_num_mons+nunkn_dep(eval_denom)+1):end] = -transpose(eval_denom.param_only[:, 1:n_instances]).*deck_values
            row_offset += n_instances
        end
    end
    return A
end

function _all_interpolated(symmetries::Vector{Vector{MiExpression}})
    for symmetry in symmetries
        for expr in symmetry
            ismissing(expr) && return false
        end
    end
    return true
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
    samples::Dict{Vector{Int}, Samples},
    path_ids_for_deck::Dict{Vector{Int}, Vector{Int}},
    eval_num_mons::Dict{Samples, EvaluatedMonomials},
    eval_denom_mons::Dict{Samples, EvaluatedMonomials},
    num_mons::AbstractMonomialVector,
    denom_mons::AbstractMonomialVector,
    tols::Tolerances;
    logging::Bool=false
)

    min_nconstraints = length(num_mons) + length(denom_mons) + 2  # TODO: x^2 + ax + b example requires (+1)
    min_ninstances = nparam_only(num_mons) + nparam_only(denom_mons)  # TODO: understand
    A = _deck_vandermonde_matrix(
        deck_permutation,
        function_id,
        samples,
        path_ids_for_deck,
        eval_num_mons,
        eval_denom_mons,
        min_nconstraints,
        min_ninstances
    )
    logging && println(
        "Created vandermonde matrix of size ",
        size(A)
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
    samples::Dict{Vector{Int}, Samples},
    path_ids_for_deck::Dict{Vector{Int}, Vector{Int}},
    eval_mons::Dict{Samples, EvaluatedMonomials},
    mons::AbstractMonomialVector,
    tols::Tolerances;
    logging::Bool=false
)

    return _interpolate_deck_function(
        deck_permutation,
        function_id,
        samples,
        path_ids_for_deck,
        eval_mons,
        eval_mons,
        mons,
        mons,
        tols;
        logging=logging
    )
end

function orbit(deck_permutations::Vector{Vector{Int}}, el)
    return unique(vcat([perm[el] for perm in deck_permutations]...))
end

function path_ids_for_deck_computation(F::SampledSystem)
    path_ids_all_deck = [Dict{Vector{Int}, Vector{Int}}() for _ in deck_permutations(F)]
    for (i, deck) in enumerate(deck_permutations(F))
        path_ids_deck = path_ids_all_deck[i]
        for path_ids_samples in keys(samples(F))
            if length(path_ids_samples) == nsolutions(F)
                path_ids_deck[path_ids_samples] = path_ids_samples
            else
                path_ids_deck[path_ids_samples] = []
                for path_id in path_ids_samples
                    if deck[path_id] in path_ids_samples
                        push!(path_ids_deck[path_ids_samples], path_id)
                    end
                end
            end
        end
    end
    return path_ids_all_deck
end

function min_nconstraints_among_deck(
    path_ids_all_deck::Vector{Dict{Vector{Int}, Vector{Int}}},
    samples::Dict{Vector{Int}, Samples}
)
    return min([sum(length(path_ids_deck)*ninstances(samples[path_ids_samples]) for (path_ids_samples, path_ids_deck) in path_ids_all_deck[i]) for i in 1:length(path_ids_all_deck)]...)
end

function _sample_for_deck_computation!(
    F::SampledSystem;
    min_nconstraints::Int,
    min_ninstances::Int
)
    path_ids_all_deck = path_ids_for_deck_computation(F)
    min_nconst = min_nconstraints_among_deck(path_ids_all_deck, samples(F))

    Δ_ninstances = min_ninstances - ninstances(F)
    Δ_min_nconst = min_nconstraints - min_nconst

    n_sols = nsolutions(F)

    if Δ_ninstances > 0
        if Δ_min_nconst > 0
            sample!(F; n_instances=div(Δ_min_nconst, n_sols))
            Δ_ninstances -= div(Δ_min_nconst, n_sols)
            if Δ_ninstances > 0
                path_ids = orbit(deck_permutations(F), 1:max(div(mod(Δ_min_nconst, n_sols), Δ_ninstances, RoundUp), 1))
                sample!(F, path_ids=path_ids, n_instances=Δ_ninstances)
            else
                sample!(F; path_ids=orbit(deck_permutations(F), 1:mod(Δ_min_nconst, n_sols)), n_instances=1)
            end
        else
            sample!(F; path_ids=orbit(deck_permutations(F), 1), n_instances=Δ_ninstances)
        end
    else
        if Δ_min_nconst > 0
            sample!(F; n_instances=div(Δ_min_nconst, n_sols))
            sample!(F; path_ids=orbit(deck_permutations(F), 1:mod(Δ_min_nconst, n_sols)), n_instances=1)
        end
    end
    return path_ids_for_deck_computation(F)
end

function symmetries_fixing_parameters_graded!(
    F::SampledSystem,
    grading::Grading;
    degree_bound::Integer=1,
    param_dep::Bool=true,
    tols::Tolerances=Tolerances(),
    logging::Bool=false
)
    C = deck_permutations(F)
    symmetries = _init_symmetries(length(C), unknowns(F))
    mons = DenseMonomialVector{Int8, Int16}(unknowns=unknowns(F), parameters=parameters(F))

    for d in 1:degree_bound
        extend!(mons; degree=d, extend_params=param_dep)
        mon_classes = to_classes(mons, grading)

        mon_classes_vect = collect(values(mon_classes))
        max_n_mons = max(length.(mon_classes_vect)...)
        max_nparam_only = max(nparam_only.(mon_classes_vect)...)

        path_ids_for_deck = _sample_for_deck_computation!(
            F;
            min_nconstraints = 3*max_n_mons+1,  # TODO: understand this
            min_ninstances = 4*max_nparam_only  # TODO: understand this
        )
        
        for (num_deg, num_mons) in mon_classes
            eval_num_mons = nothing
            for i in 1:nunknowns(F)
                denom_deg = _denom_deg(num_deg, grading, i)  # i-th variable
                denom_mons = get(mon_classes, denom_deg, nothing)
                if !isnothing(denom_mons)
                    num_denom_mons = vcat(num_mons, denom_mons)
                    if iszero(gcd(num_denom_mons)) && !is_param_only(num_denom_mons)
                        if isnothing(eval_num_mons)
                            eval_num_mons = evaluate(num_mons, samples(F))
                        end
                        eval_denom_mons = evaluate(denom_mons, samples(F))
                        for (j, symmetry) in enumerate(symmetries)
                            if ismissing(symmetry[i])
                                symmetry[i] = _interpolate_deck_function(
                                    C[j],
                                    i,
                                    samples(F),
                                    path_ids_for_deck[j],
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
    end

    return DeckTransformationGroup(symmetries, F)
end

function symmetries_fixing_parameters_dense!(
    F::SampledSystem; 
    degree_bound::Integer=1,
    param_dep::Bool=true,
    tols::Tolerances=Tolerances(),
    logging::Bool=false
)
    C = deck_permutations(F)
    symmetries = _init_symmetries(length(C), unknowns(F))
    mons = DenseMonomialVector{Int8, Int16}(unknowns=unknowns(F), parameters=parameters(F))

    for d in 1:degree_bound
        logging && printstyled("Started interpolation for degree = ", d, "...\n"; color=:green)
        
        extend!(mons; degree=d, extend_params=param_dep)
        path_ids_for_deck = _sample_for_deck_computation!(
            F;
            min_nconstraints = 2*length(mons)+1,  # TODO: understand this
            min_ninstances = 2*nparam_only(mons)  # TODO: understand this
        )

        logging && println("Evaluating monomials...\n")
        # TODO: pick samples for evaluation
        evaluated_mons = evaluate(mons, samples(F))

        for (i, symmetry) in enumerate(symmetries)
            logging && printstyled("Interpolating the ", i, "-th deck transformation...\n"; color=:blue)
            for j in 1:nunknowns(F)
                if ismissing(symmetry[j])
                    symmetry[j] = _interpolate_deck_function(
                        C[i],
                        j,
                        samples(F),
                        path_ids_for_deck[i],
                        evaluated_mons,
                        mons,
                        tols;
                        logging=logging
                    )
                    if logging && !ismissing(symmetry[j])
                        printstyled(
                            "Good representative for the ",
                            i,
                            "-th deck transformation, variable ",
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
            logging && printstyled("--- All deck transformations are interpolated ---\n"; color=:blue)
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
    final_grading = Grading(nfree(grading), grading.free_part, [])
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
    symmetries_fixing_parameters(F::System; degree_bound=1, param_dep=true, kwargs...)

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
   x₁ ↦ -x₁
   x₂ ↦ -x₂
  3rd map:
   x₁ ↦ im*x₂
   x₂ ↦ -im*x₁
  4th map:
   x₁ ↦ -im*x₂
   x₂ ↦ im*x₁
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

# ----------------------- DEBUGGING TOOLS -----------------------
export deck_vandermonde_dense,
    deck_vandermonde_graded,
    to_multiexponent,
    reduced_nullspace,
    to_coefficients

function reduced_nullspace(A::AbstractMatrix{<:Number}; tols::Tolerances=Tolerances())
    if tols.nullspace_rtol == 0
        N = nullspace(A, atol=tols.nullspace_atol)
    else
        N = nullspace(A, atol=tols.nullspace_atol, rtol=tols.nullspace_rtol)
    end
    N = transpose(N)
    size(N, 1) == 0 && return N
    N = rref(N, tols.rref_tol)
    sparsify!(N, tols.sparsify_tol; digits=1)
    return N
end

function deck_vandermonde_dense(
    F::SampledSystem;
    deck_id::Int,
    var::Variable,
    degree::Integer=1,
    param_dep::Bool=true,
    min_nconstraints::Function=mons->2*length(mons),
    min_ninstances::Function=mons->2*nparam_only(mons)
)
    mons = DenseMonomialVector{Int8, Int16}(unknowns=unknowns(F), parameters=parameters(F))
    extend!(mons; degree=degree, extend_params=param_dep)
    path_ids_for_deck = _sample_for_deck_computation!(
        F;
        min_nconstraints = min_nconstraints(mons),
        min_ninstances = min_ninstances(mons)
    )
    eval_mons = evaluate(mons, samples(F))
    function_id = findfirst(x->x==var, unknowns(F))
    return _deck_vandermonde_matrix(
        deck_permutations(F)[deck_id],
        function_id,
        samples(F),
        path_ids_for_deck[deck_id],
        eval_mons,
        eval_mons,
        min_nconstraints(mons),
        min_ninstances(mons)
    ), mons
end

function to_multiexponent(mon::Expression, vars::Vector{Variable})
    es, cs = exponents_coefficients(mon, vars)
    if length(cs) > 1 || !isone(cs[1])
        throw(ArgumentError("Input expression is not a monomial"))
    else
        return SparseVector{Int8,Int16}(sparse(es[:,1]))
    end
end

function deck_vandermonde_graded(
    F::SampledSystem,
    grading::Grading;
    deck_id::Int,
    var::Variable,
    deg::Integer=1,
    param_dep::Bool=true,
    num_mon::Expression,
    min_nconstraints::Function=(num_mons, denom_mons)->length(num_mons)+length(denom_mons)+2,
    min_ninstances::Function=(num_mons, denom_mons)->nparam_only(num_mons)+nparam_only(denom_mons)
)
    mons = DenseMonomialVector{Int8, Int16}(unknowns=unknowns(F), parameters=parameters(F))
    extend!(mons; degree=deg, extend_params=param_dep)
    mon_classes = to_classes(mons, grading)

    var_id = findfirst(x->x==var, variables(F))

    mexp = to_multiexponent(num_mon, variables(F))
    num_deg = degree(mexp, grading)
    num_mons = mon_classes[num_deg]
    denom_deg = _denom_deg(num_deg, grading, var_id)
    denom_mons = get(mon_classes, denom_deg, nothing)
    isnothing(denom_mons) && error("denom_mons don't exist")

    path_ids_for_deck = _sample_for_deck_computation!(
        F;
        min_nconstraints = min_nconstraints(num_mons, denom_mons),
        min_ninstances = min_ninstances(num_mons, denom_mons)
    )

    eval_num_mons = evaluate(num_mons, samples(F))
    eval_denom_mons = evaluate(denom_mons, samples(F))
    return _deck_vandermonde_matrix(
        deck_permutations(F)[deck_id],
        var_id,
        samples(F),
        path_ids_for_deck[deck_id],
        eval_num_mons,
        eval_denom_mons,
        min_nconstraints(num_mons, denom_mons),
        min_ninstances(num_mons, denom_mons)
    ), (num_mons, denom_mons)
end

function to_coefficients(
    expr::Expression,
    mons::SparseMonomialVector{Tv,Ti}
) where {Tv<:Integer,Ti<:Integer}
    es, cs = exponents_coefficients(expr, variables(mons))
    mexps = [SparseVector{Tv,Ti}(sparse(esᵢ)) for esᵢ in eachcol(es)]
    ids = [findfirst(m, mons) for m in mexps]
    v = zeros(ComplexF64, length(mons))
    v[ids] = cs
    return v
end

function to_coefficients(
    num_expr::Expression,
    denom_expr::Expression,
    num_mons::SparseMonomialVector{Tv,Ti},
    denom_mons::SparseMonomialVector{Tv,Ti}
) where {Tv<:Integer,Ti<:Integer}
    return vcat(to_coefficients(num_expr, num_mons), to_coefficients(denom_expr, denom_mons))
end

function to_coefficients(
    num_expr::Expression,
    denom_expr::Expression,
    mons::SparseMonomialVector{Tv,Ti}
) where {Tv<:Integer,Ti<:Integer}
    return to_coefficients(num_expr, denom_expr, mons, mons)
end