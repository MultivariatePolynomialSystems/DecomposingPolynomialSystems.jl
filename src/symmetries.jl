include("interpolation.jl")

export DeckTransformationGroup, ScalingSymmetryGroup
export show_map
export scaling_symmetries, _verify_commutativity
export symmetries_fixing_parameters_dense!, symmetries_fixing_parameters_graded!
export symmetries_fixing_parameters!, symmetries_fixing_parameters

struct DeckTransformationGroup
    maps::Vector{RationalMap}
    structure::String
    F::SampledSystem
end

function DeckTransformationGroup(F::SampledSystem)
    symmetries = _init_symmetries(length(F.deck_permutations), unknowns(F))
    return DeckTransformationGroup(symmetries, F)
end

function DeckTransformationGroup(
    symmetries::Vector{Vector{NoExpression}},
    F::SampledSystem
)
    action = [RationalMap(variables(F), vcat(symmetry, parameters(F))) for symmetry in symmetries]
    return DeckTransformationGroup(action, group_structure(F.deck_permutations), F)
end

function Base.show(io::IO, deck::DeckTransformationGroup)
    println(io, "DeckTransformationGroup of order $(length(deck.maps))")
    println(io, " structure: ", deck.structure)
    print(io, " action:")
    for i in eachindex(deck.maps)
        println(io, "\n  ", to_ordinal(i), " map:")
        for (j, var) in enumerate(unknowns(deck.F))  # action on parameters is trivial, don't show it
            print(io, "   ", var, " ↦ ", deck.maps[i][var])
            j < length(unknowns(deck.F)) && print(io, "\n")
        end
    end
end

Base.getindex(deck::DeckTransformationGroup, i::Int) = deck.maps[i]

# TODO: make pretty printing
struct ScalingSymmetryGroup
    grading::Grading
    vars::Vector{Variable}
    action::Vector{Tuple{Int, Dict{Variable, Expression}}}
end

ScalingSymmetryGroup() = ScalingSymmetryGroup([], [], [])

function ScalingSymmetryGroup(grading::Grading, vars::Vector{Variable})
    action = Vector{Tuple{Int, Dict{Variable, Expression}}}([])
    k = 0
    for (sᵢ, Uᵢ) in grading
        n_scalings = size(Uᵢ, 1)
        @var λ[(k+1):(n_scalings+k)]
        for j in 1:n_scalings
            nonzero_ids = findall(!iszero, Uᵢ[j, :])
            vals = (λ[j].^Uᵢ[j, :][nonzero_ids]).*vars[nonzero_ids]
            push!(action, (sᵢ, Dict(zip(vars[nonzero_ids], vals))))
        end
        k += n_scalings
    end
    return ScalingSymmetryGroup(grading, vars, action)
end

Base.copy(s::ScalingSymmetryGroup) = ScalingSymmetryGroup(s.grading, s.vars, s.action)

function Base.show(io::IO, scalings::ScalingSymmetryGroup)
    n_infinite, n_finite = 0, 0
    for (sᵢ, Uᵢ) in scalings.grading
        if sᵢ == 0
            n_infinite = size(Uᵢ, 1)
        else
            n_finite += size(Uᵢ, 1)
        end
    end
    println(io, "ScalingSymmetryGroup with $(n_infinite+n_finite) scalings")
    println(io, " infinite scalings: $(n_infinite)")
    println(io, " finite scalings:")
    for (i, (sᵢ, Uᵢ)) in enumerate(scalings.grading)
        if sᵢ == 0 continue end
        println(io, "  $(size(Uᵢ, 1)) of order $(sᵢ)")
    end
    print(io, " vars: ", join(scalings.vars, ", "))
end

function _mat2col_diffs(M)
    M = M - M[:,1]*ones(eltype(M), 1, size(M,2))
    return M[:,2:end]
end

function _snf_scaling_symmetries(F::System)::Tuple{Vector{Int}, Vector{Matrix{Int}}}
    vars = vcat(F.variables, F.parameters)
    Es = [exponents_coefficients(f, vars)[1] for f in F.expressions]
    K = VM2M([_mat2col_diffs(E) for E in Es])
    if size(K, 1) > size(K, 2)
        K = [K zeros(eltype(K), size(K, 1), size(K,1)-size(K,2))]
    end
    K = matrix(ZZ, K)

    S, U, _ = snf_with_transform(K)
    U, S = Matrix(U), Int.(diag(Matrix(S)))

    s = reverse(filter(el->el!=1, unique(S)))
    if length(s) == 0
        return [], []
    end
    idxs = [findall(x->x==el, S) for el in s]
    Us = Vector{Matrix{Int}}([])
    for i in eachindex(idxs)
        push!(Us, U[idxs[i], :])
    end
    return s, Us # TODO: what if max(Us[i]) > MAX_INT64?
end

function _hnf_reduce!(grading::Grading)
    for (i, (sᵢ, Uᵢ)) in enumerate(grading)
        if sᵢ == 0
            grading[i] = (sᵢ, Matrix(hnf(matrix(ZZ, Uᵢ))))
        else
            grading[i] = (sᵢ, lift.(Matrix(hnf(matrix(GF(sᵢ), Uᵢ)))))
        end
    end
end

"""
    scaling_symmetries(F::System; in_hnf::Bool=true)

Given a polynomial system `F` returns the group of scaling symmetries 
of the polynomial system `F`.

```julia-repl
julia> @var x[1:2] p[1:2];

julia> F = System([x[1]^2 - x[2]^2 - p[1], 2*x[1]*x[2] - p[2]]; variables=x, parameters=p);

julia> scalings = scaling_symmetries(F)
ScalingSymmetryGroup with 2 scalings
 infinite scalings: 1
 finite scalings:
  1 of order 2
 vars: x₁, x₂, p₁, p₂
```
"""
function scaling_symmetries(F::System; in_hnf::Bool=true)::ScalingSymmetryGroup
    s, U = _snf_scaling_symmetries(F)
    if length(s) == 0
        return ScalingSymmetryGroup()
    end
    grading = collect(zip(s, U))
    if in_hnf _hnf_reduce!(grading) end
    vars = vcat(F.variables, F.parameters)
    return ScalingSymmetryGroup(grading, vars)
end

function scaling_symmetries(F::SampledSystem; in_hnf::Bool=true)::ScalingSymmetryGroup
    return scaling_symmetries(F.system; in_hnf=in_hnf)
end

function _remove_zero_rows(M::Matrix)::Matrix
    nonzero_rows = filter(!iszero, M2VV(transpose(M)))
    if length(nonzero_rows) == 0
        return Matrix{eltype(M)}(undef, 0, size(M, 2))
    end
    return transpose(VV2M(nonzero_rows))
end

function _remove_zero_rows(grading::Grading)::Grading
    filtered_grading = Grading([])
    for (i, (sᵢ, Uᵢ)) in enumerate(grading)
        Uᵢ = _remove_zero_rows(Uᵢ)
        if size(Uᵢ, 1) != 0
            push!(filtered_grading, (sᵢ, Uᵢ))
        end
    end
    return filtered_grading
end

function _remove_dependencies(grading::Grading)::Grading
    # TODO: remove rows dependent on other blocks
    return grading
end

# TODO: what if vars is not a subset of scalings.vars?
function _restrict_scalings(scalings::ScalingSymmetryGroup, vars::Vector{Variable})::ScalingSymmetryGroup
    idx = [findfirst(v->v==var, scalings.vars) for var in vars]
    restr_grading = copy(scalings.grading)
    for (i, (sᵢ, Uᵢ)) in enumerate(scalings.grading)
        restr_grading[i] = (sᵢ, Uᵢ[:, idx])
    end
    _hnf_reduce!(restr_grading)
    restr_grading = _remove_dependencies(_remove_zero_rows(restr_grading))
    return ScalingSymmetryGroup(restr_grading, vars)
end

function scaling_symmetries(F::System, vars::Vector{Variable})::ScalingSymmetryGroup
    return _restrict_scalings(scaling_symmetries(F), vars)
end

function scaling_symmetries(F::SampledSystem, vars::Vector{Variable})::ScalingSymmetryGroup
    return scaling_symmetries(F.system, vars)
end

function _num_deg2denom_deg(num_deg::Vector{Int}, grading::Grading, var_id::Int)::Vector{Int}
    denom_deg = zeros(Int, length(num_deg))
    k = 1
    for (sᵢ, Uᵢ) in grading
        n_scalings = size(Uᵢ, 1)
        denom_deg[k:k+n_scalings-1] = modV(num_deg[k:k+n_scalings-1] - Uᵢ[:, var_id], sᵢ)
        k += n_scalings
    end
    return denom_deg
end

# TODO: Can it be eps close to zero? Then the method isn't correctly written...
# TODO: Do we need printing? Then passing mons isn't necessary, just their number
# TODO: change the name of the method?
function _remove_zero_nums_and_denoms(
    coeffs::Matrix{CC},
    num_mons::MonomialVector,
    denom_mons::MonomialVector;
    logging::Bool=false
)::Matrix{CC}

    reasonable_rows = []
    n_num_mons, n_denom_mons = length(num_mons.mds), length(denom_mons.mds)
    @assert size(coeffs, 2) == n_num_mons + n_denom_mons
    for i in 1:size(coeffs, 1)
        if (!all(iszero, coeffs[i, 1:n_num_mons]) && !all(iszero, coeffs[i, n_num_mons+1:end]))
            push!(reasonable_rows, i)
        elseif logging
            println("Removed: ",
                dot(coeffs[i, 1:n_num_mons], num_mons) + dot(coeffs[i, n_num_mons+1:end], denom_mons)  # TODO: convert to MonomialVector input
            )
        end
    end
    return coeffs[reasonable_rows, :]
end

function _remove_zero_nums_and_denoms(
    coeffs::Matrix{CC},
    mons::MonomialVector;
    logging::Bool=false
)::Matrix{CC}

    return _remove_zero_nums_and_denoms(coeffs, mons, mons, logging=logging)
end

function _vandermonde_matrix(
    permutation::Vector{Int},
    values::SubArray{CC, 2},
    eval_num_mons::Array{CC, 3},
    eval_denom_mons::Array{CC, 3}
)::Matrix{CC}

    n_sols, n_instances = size(values)
    n_num_mons = size(eval_num_mons, 1)
    n_denom_mons = size(eval_denom_mons, 1)

    A = zeros(CC, n_instances*n_sols, n_num_mons+n_denom_mons)
    @assert size(A, 1) >= size(A, 2)

    for i in 1:n_sols
        v = values[permutation[i], :]
        rows = ((i-1)*n_instances+1):(i*n_instances)
        A[rows, 1:n_num_mons] = transpose(eval_num_mons[:, i, :])
        A[rows, (n_num_mons+1):end] = -transpose(eval_denom_mons[:, i, :]).*v
    end
    return A
end

function _vandermonde_matrix(
    permutation::Vector{Int}, 
    values::SubArray{CC, 2}, 
    eval_mons::Array{CC, 3}
)::Matrix{CC}

    return _vandermonde_matrix(permutation, values, eval_mons, eval_mons)
end

function _all_interpolated(symmetries::Vector{Vector{NoExpression}})::Bool
    all_interpolated = true
    for symmetry in symmetries
        if nothing in symmetry
            all_interpolated = false
            break
        end
    end
    return all_interpolated
end

function _init_symmetries(n_symmetries::Int, unknowns::Vector{Variable})::Vector{Vector{NoExpression}}
    symmetries = [[nothing for j in eachindex(unknowns)] for i in 1:n_symmetries]
    symmetries = Vector{Vector{NoExpression}}(symmetries)
    symmetries[1] = Expression.(unknowns)  # set the first to the identity
    return symmetries
end

function _interpolate_symmetry_function(
    permutation::Vector{Int},
    values::SubArray{CC, 2},
    eval_num_mons::Array{CC, 3},
    eval_denom_mons::Array{CC, 3},
    num_mons::MonomialVector,
    denom_mons::MonomialVector,
    tol::Float64;
    logging::Bool=false
)::NoExpression

    logging && println(
        "Creating vandermonde matrix of size ",
        (prod(size(values)), length(num_mons)+length(denom_mons))
    )
    A = _vandermonde_matrix(permutation, values, eval_num_mons, eval_denom_mons)

    logging && println("Computing nullspace...")
    coeffs = Matrix{CC}(transpose(nullspace(A)))
    logging && println("Size of the transposed nullspace: ", size(coeffs))

    if size(coeffs, 1) == 0
        return nothing
    end

    logging && println("Computing the reduced row echelon form of the transposed nullspace...\n")
    coeffs = sparsify(rref(coeffs, tol), tol, digits=1)
    coeffs = _remove_zero_nums_and_denoms(coeffs, num_mons, denom_mons)

    if size(coeffs, 1) == 0
        return nothing
    end
    coeffs = good_representative(coeffs)
    return rational_function(coeffs, num_mons, denom_mons; logging=false, tol=tol)
end

function _interpolate_symmetry_function(
    permutation::Vector{Int},
    values::SubArray{CC, 2},
    eval_mons::Array{CC, 3},
    mons::MonomialVector,
    tol::Float64;
    logging::Bool=false
)

    return _interpolate_symmetry_function(
        permutation,
        values,
        eval_mons,
        eval_mons,
        mons,
        mons,
        tol;
        logging=logging
    )
end

function symmetries_fixing_parameters_graded!(
    F::SampledSystem,
    scalings::ScalingSymmetryGroup,
    mds::Vector{Multidegree},
    classes::Dict{Vector{Int}, Vector{Int}};
    tol::Float64=1e-5,
    logging::Bool=false
)::DeckTransformationGroup
    
    max_n_mons = max(length.(collect(values(classes)))...)  # size of the largest class
    n_unknowns, n_sols, _ = size(F.samples.solutions)  # TODO: what if n_sols is huge?
    n_instances = Int(ceil(2/n_sols*max_n_mons))

    C = F.deck_permutations
    symmetries = _init_symmetries(length(C), unknowns(F))

    sample_system!(F, n_instances)
    
    for (num_deg, num_ids) in classes
        num_mons = MonomialVector(mds[num_ids], scalings.vars)
        eval_num_mons = nothing
        for i in 1:n_unknowns
            denom_deg = _num_deg2denom_deg(num_deg, scalings.grading, i)  # i-th variable
            denom_ids = get(classes, denom_deg, nothing)
            if !isnothing(denom_ids)
                denom_mons = MonomialVector(mds[denom_ids], scalings.vars)
                g = gcd_mds([gcd_mons(num_mons), gcd_mons(denom_mons)])  # TODO: improve
                if iszero(g) && (!only_param_dep(num_mons, n_unknowns) || !only_param_dep(denom_mons, n_unknowns))  # TODO: improve
                    if isnothing(eval_num_mons)
                        eval_num_mons = evaluate_monomials_at_samples_(num_mons, F.samples)
                    end
                    eval_denom_mons = evaluate_monomials_at_samples_(denom_mons, F.samples)
                    for (j, symmetry) in enumerate(symmetries)
                        if isnothing(symmetry[i])
                            symmetry[i] = _interpolate_symmetry_function(
                                C[j],
                                view(F.samples.solutions, i, :, :),
                                eval_num_mons,
                                eval_denom_mons,
                                num_mons,
                                denom_mons,
                                tol;
                                logging=logging
                            )
                            if !isnothing(symmetry[i])
                                logging && printstyled(
                                    "Good representative for the ",
                                    to_ordinal(j),
                                    " symmetry, variable ",
                                    unknowns(F)[i],
                                    ":\n",
                                    color=:red
                                )
                                logging && println(symmetry[i])
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
    scalings::ScalingSymmetryGroup;
    degree_bound::Int=1,
    tol::Float64=1e-5,
    logging::Bool=false
)::DeckTransformationGroup

    mds = multidegrees_up_to_total_degree(length(scalings.vars), degree_bound)
    classes = partition_multidegrees(mds, scalings.grading)
    return symmetries_fixing_parameters_graded!(
        F,
        scalings,
        mds,
        classes;
        tol=tol,
        logging=logging
    )
end

function symmetries_fixing_parameters_dense!(
    F::SampledSystem; 
    degree_bound::Int=1,
    param_dep::Bool=true,
    tol::Float64=1e-5,
    logging::Bool=false
)::DeckTransformationGroup

    n_unknowns, n_sols, _ = size(F.samples.solutions)  # TODO: what if n_sols is huge?
    vars = param_dep ? variables(F) : unknowns(F)  # rename vars --> interp_vars?

    C = F.deck_permutations
    symmetries = _init_symmetries(length(C), unknowns(F))

    for d in 1:degree_bound
        logging && printstyled("Started interpolation for degree = ", d, "...\n", color=:green)
        mons = monomials(vars, d)
        n_instances = Int(ceil(2/n_sols*length(mons)))
        sample_system!(F, n_instances)

        logging && println("Evaluating monomials...\n")
        evaluated_mons = evaluate_monomials_at_samples_(mons, F.samples)
        
        for (i, symmetry) in enumerate(symmetries)
            logging && printstyled("Interpolating the ", i, "-th symmetry map...\n", color=:blue)
            for j in 1:n_unknowns
                if isnothing(symmetry[j])
                    symmetry[j] = _interpolate_symmetry_function(
                        C[i],
                        view(F.samples.solutions, j, :, :),
                        evaluated_mons,
                        mons,
                        tol;
                        logging=logging
                    )
                    if !isnothing(symmetry[j])
                        logging && printstyled(
                            "Good representative for the ",
                            i,
                            "-th symmetry, variable ",
                            unknowns(F)[j],
                            ":\n",
                            color=:red
                        )
                        logging && println(symmetry[j])
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

function to_CC(scaling::Tuple{Int, Vector{Int}})::Vector{CC}
    return [exp(2*pi*im*k/scaling[1]) for k in scaling[2]]
end

function find_solution_id(sol::Vector{CC}, sols::Matrix{CC}, tol::Float64)
    return findfirst(x->norm(x-sol)<tol, M2VV(sols))
end

# verify for all of the solutions in 1 instance
function _all_deck_commute(F::SampledSystem, scaling::Tuple{Int, Vector{Int}}; tol::Float64=1e-5)::Bool
    instance_id = rand(1:size(F.samples.solutions, 3))
    sols1 = F.samples.solutions[:, :, instance_id]
    params1 = F.samples.parameters[:, instance_id]
    params2 = to_CC(scaling)[end-n_parameters(F)+1:end].*params1
    println("Starting sampling...")
    sample_system!(F, params2)  # TODO: what if n_sols is huge?
    println("Finished sampling...")
    sols2 = F.samples.solutions[:, :, end]
    for perm in F.deck_permutations
        for i in axes(sols1, 2)
            Φ_sol = to_CC(scaling)[1:n_unknowns(F)].*sols1[:, i]
            id = find_solution_id(Φ_sol, sols2, tol)
            if isnothing(id) return false end
            ΨΦ_sol = sols2[:, perm[id]]
            Ψ_sol = sols1[:, perm[i]]
            ΦΨ_sol = to_CC(scaling)[1:n_unknowns(F)].*Ψ_sol
            if norm(ΨΦ_sol-ΦΨ_sol)>tol return false end
        end
    end
    return true
end

# TODO: consider every element from finite components => take every linear combination with coeffs from Z_si
function _scalings_commuting_with_deck(F::SampledSystem, scalings::ScalingSymmetryGroup)
    final_grading = Grading([])
    for (sᵢ, Uᵢ) in scalings.grading
        if sᵢ == 0
            push!(final_grading, (sᵢ, Uᵢ))
            continue
        end
        Vᵢ = Array{Int}(undef, 0, size(Uᵢ, 2))
        for j in axes(Uᵢ, 1)
            if _all_deck_commute(F, (sᵢ, Uᵢ[j, :]))
                Vᵢ = [Vᵢ; V2Mt(Uᵢ[j, :])]
            end
        end
        if size(Vᵢ, 1) > 0
            push!(final_grading, (sᵢ, Vᵢ))
        end
    end
    return ScalingSymmetryGroup(final_grading, scalings.vars)
end

function symmetries_fixing_parameters!(
    F::SampledSystem;
    degree_bound::Int=1,
    param_dep::Bool=true,
    tol::Float64=1e-5,
    logging::Bool=false
)::DeckTransformationGroup

    if length(F.deck_permutations) == 1 # trivial group of symmetries
        return DeckTransformationGroup(F) # return the identity group
    end

    # scalings = scaling_symmetries(F)
    scalings = _scalings_commuting_with_deck(F, scaling_symmetries(F))
    scalings = param_dep ? scalings : _restrict_scalings(scalings, unknowns(F))
    if length(scalings.grading) == 0
        return symmetries_fixing_parameters_dense!(
            F;
            degree_bound=degree_bound,
            param_dep=param_dep,
            tol=tol,
            logging=logging
        )
    else
        logging && printstyled("Running graded version...\n", color=:green)
        return symmetries_fixing_parameters_graded!(
            F,
            scalings;
            degree_bound=degree_bound,
            tol=tol,
            logging=logging
        )
    end
end

function symmetries_fixing_parameters(  # TODO: extend to take a rational map
    F::System,
    (x₀, p₀)::Tuple{Vector{CC}, Vector{CC}};
    degree_bound::Int=1,
    param_dep::Bool=true,
    tol::Float64=1e-5,
    monodromy_options::Tuple=(),
    logging::Bool=false
)::DeckTransformationGroup

    F = run_monodromy(F, (x₀, p₀); monodromy_options...)
    return symmetries_fixing_parameters!(
        F;
        degree_bound=degree_bound,
        param_dep=param_dep,
        tol=tol,
        logging=logging
    )
end

"""
    symmetries_fixing_parameters(F::System; degree_bound=1, param_dep=true, tol=1e-5)

Given a polynomial system F returns the group of symmetries 
of the polynomial system `F` that fix the parameters. The keyword
argument `degree_bound` is used to set the upper bound for the
degrees of numerator and denominator polynomials in expressions
for the symmetries.

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
function symmetries_fixing_parameters(  # TODO: extend to take a rational map
    F::System;
    degree_bound::Int=1,
    param_dep::Bool=true,
    tol::Float64=1e-5,
    monodromy_options::Tuple=(),
    logging::Bool=false
)::DeckTransformationGroup

    F = run_monodromy(F; monodromy_options...)
    return symmetries_fixing_parameters!(
        F;
        degree_bound=degree_bound,
        param_dep=param_dep,
        tol=tol,
        logging=logging
    )
end

