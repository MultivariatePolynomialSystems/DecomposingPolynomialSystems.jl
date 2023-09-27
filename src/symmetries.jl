include("interpolation.jl")

export DeckTransformationGroup, ScalingSymmetryGroup
export scaling_symmetries
export symmetries_fixing_parameters_dense!, symmetries_fixing_parameters_graded!
export symmetries_fixing_parameters!, symmetries_fixing_parameters

NoExpression = Union{Nothing, Expression}

struct RationalMap

end

struct DeckTransformationGroup
    action::Vector{Dict{Variable, NoExpression}}
    F::SampledSystem
    # structure::GapObj
end

function DeckTransformationGroup(F::SampledSystem)
    action = [Dict(zip(unknowns(F), Expression.(unknowns(F))))]  # TODO: Is conversion needed?
    return DeckTransformationGroup(action, F)
end

function DeckTransformationGroup(
    symmetries::Vector{Vector{NoExpression}},
    F::SampledSystem
)
    action = [Dict(zip(unknowns(F), symmetry)) for symmetry in symmetries]
    return DeckTransformationGroup(action, F)
end

struct ScalingSymmetryGroup
    grading::Grading
    vars::Vector{Variable}
    action::Vector{Tuple{Int, Dict{Variable, Expression}}}
end

ScalingSymmetryGroup() = ScalingSymmetryGroup([], [], [])
Base.copy(s::ScalingSymmetryGroup) = ScalingSymmetryGroup(s.grading, s.vars, s.action)

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

function _to_expressions(grading::Grading, vars::Vector{Variable})
    action = Vector{Tuple{Int, Dict{Variable, Expression}}}([])
    k = 0
    for (sᵢ, Uᵢ) in grading
        n_scalings = size(Uᵢ, 1)
        @var λ[(k+1):(n_scalings+k)]
        for j in 1:n_scalings
            nonzero_ids = findall(!iszero, Uᵢ[j, :])
            println("k = ", k)
            println("j = ", j)
            vals = (λ[j].^Uᵢ[j, :][nonzero_ids]).*vars[nonzero_ids]
            push!(action, (sᵢ, Dict(zip(vars[nonzero_ids], vals))))
        end
        k += n_scalings
    end
    return action
end

function scaling_symmetries(F::System; in_hnf::Bool=true)::ScalingSymmetryGroup
    s, U = _snf_scaling_symmetries(F)
    if length(s) == 0
        return ScalingSymmetryGroup()
    end
    grading = collect(zip(s, U))
    if in_hnf _hnf_reduce!(grading) end
    vars = vcat(F.variables, F.parameters)
    action = _to_expressions(grading, vars)
    return ScalingSymmetryGroup(grading, vars, action)
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

# TODO: what if vars is not a subset of Fvars?
function scaling_symmetries(F::System, vars::Vector{Variable})::ScalingSymmetryGroup
    Fvars = vcat(F.variables, F.parameters)
    idx = [findfirst(v->v==var, Fvars) for var in vars]
    s, U = _snf_scaling_symmetries(F)
    if length(s) == 0
        return ScalingSymmetryGroup()
    end
    U = [u[:, idx] for u in U]
    grading = collect(zip(s, U))
    _hnf_reduce!(grading)
    grading = _remove_dependencies(_remove_zero_rows(grading))
    action = _to_expressions(grading, vars)
    return ScalingSymmetryGroup(grading, vars, action)
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
    if all_interpolated
        return true
    end
    return false
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
    mons::Vector{Expression},
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
)::DeckTransformationGroup
    
    max_n_mons = max(length.(collect(values(classes)))...)  # size of the largest class
    n_unknowns, n_sols, _ = size(F.samples.solutions)  # TODO: what if n_sols is huge?
    n_instances = Int(ceil(2/n_sols*max_n_mons))

    C = F.symmetry_permutations
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
                                tol
                            )
                            if !isnothing(symmetry[i])
                                printstyled(
                                    "Good representative for the ",
                                    j,
                                    "-th symmetry, variable ",
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
            printstyled("--- All symmetries are interpolated ---\n", color=:blue)
            # printstyled("Number of processed classes: ", m, " out of ", length(classes), "\n", color=:blue)
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
)::DeckTransformationGroup

    mds = multidegrees_up_to_total_degree(length(scalings.vars), degree_bound)
    classes = partition_multidegrees(mds, scalings.grading)
    return symmetries_fixing_parameters_graded!(
        F,
        scalings,
        mds,
        classes;
        tol=tol
    )
end

function symmetries_fixing_parameters_dense!(
    F::SampledSystem; 
    degree_bound::Int=1,
    param_dep::Bool=true,
    tol::Float64=1e-5
)::DeckTransformationGroup

    n_unknowns, n_sols, _ = size(F.samples.solutions)  # TODO: what if n_sols is huge?
    vars = param_dep ? variables(F) : unknowns(F)  # vars --> interp_vars?
    n_vars = length(vars)

    C = F.symmetry_permutations
    symmetries = _init_symmetries(length(C), unknowns(F))

    mons = Expression.([1])

    for d in 1:degree_bound
        printstyled("Started interpolation for degree = ", d, "...\n", color=:green)

        new_n_instances = Int(ceil(2/n_sols*num_mons_upto(n_vars, d)))
        new_mons = next_deg_mons(vars, mons, d)
        
        sample_system!(F, new_n_instances)
        # append!(mons, new_mons)
        prepend!(mons, new_mons)

        println("Evaluating new monomials...\n")
        evaluated_mons = evaluate_monomials_at_samples(mons, F.samples, vars) # TODO: optimize
        
        for i in 2:length(C)  # skip the identity permutation
            printstyled("Interpolating the ", i, "-th symmetry map...\n", color=:blue)
            symmetry = symmetries[i]
            for j in 1:n_unknowns
                if isnothing(symmetry[j])
                    symmetry[j] = _interpolate_symmetry_function(
                        C[i],
                        view(F.samples.solutions, j, :, :),
                        evaluated_mons,
                        mons,
                        tol;
                        logging=true
                    )
                    if !isnothing(symmetry[j])
                        printstyled(
                            "Good representative for the ",
                            i,
                            "-th symmetry, variable ",
                            unknowns(F)[j],
                            ":\n",
                            color=:red
                        )
                        println(symmetry[j])
                    end
                end
            end
        end
    
        if _all_interpolated(symmetries)
            printstyled("--- All symmetries are interpolated ---\n", color=:blue)
            return DeckTransformationGroup(symmetries, F)
        end
    end

    return DeckTransformationGroup(symmetries, F)
end

function symmetries_fixing_parameters!(
    F::SampledSystem;
    degree_bound::Int=1,
    param_dep::Bool=true,
    tol::Float64=1e-5
)::DeckTransformationGroup

    if length(F.symmetry_permutations) == 1 # trivial group of symmetries
        return DeckTransformationGroup(F) # return the identity group
    end

    scalings = param_dep ? scaling_symmetries(F) : scaling_symmetries(F, unknowns(F))
    # TODO: Verify, if finite scalings commute
    if length(scalings.grading) == 0
        return symmetries_fixing_parameters_dense!(
            F;
            degree_bound=degree_bound,
            param_dep=param_dep,
            tol=tol
        )
    else
        println("Found non-trivial grading:")
        for (sᵢ, Uᵢ) in scalings.grading
            printstyled(size(Uᵢ, 1), " scalings of order ", sᵢ, "\n", color=:green)
        end
        printstyled("Running graded version...\n", color=:green)
        return symmetries_fixing_parameters_graded!(
            F,
            scalings;
            degree_bound=degree_bound,
            tol=tol
        )
    end
end

function symmetries_fixing_parameters(  # TODO: extend to take unknowns arg
    F::System,
    (x₀, p₀)::Tuple{Vector{CC}, Vector{CC}};
    degree_bound::Int=1,
    param_dep::Bool=true,
    tol::Float64=1e-5,
    monodromy_options::Tuple=()
)::DeckTransformationGroup

    F = run_monodromy(F, (x₀, p₀); monodromy_options...)
    return symmetries_fixing_parameters!(
        F;
        degree_bound=degree_bound,
        param_dep=param_dep,
        tol=tol
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
julia> @var x[1:2] p[1:2]
(Variable[x₁, x₂], Variable[p₁, p₂])

julia> F = System([x[1]^2 - x[2]^2 - p[1], 2*x[1]*x[2] - p[2]]; variables=x, parameters=p)
System of length 2
 2 variables: x₁, x₂
 2 parameters: p₁, p₂

 -p₁ + x₁^2 - x₂^2
 -p₂ + 2*x₂*x₁

julia> deck = symmetries_fixing_parameters(F, degree_bound=1, param_dep=false)

```
"""
function symmetries_fixing_parameters(  # TODO: extend to take a rational map
    F::System;
    degree_bound::Int=1,
    tol::Float64=1e-5,
    param_dep::Bool=true,
    monodromy_options::Tuple=()
)::DeckTransformationGroup

    F = run_monodromy(F; monodromy_options...)
    return symmetries_fixing_parameters!(
        F;
        degree_bound=degree_bound,
        tol=tol,
        param_dep=param_dep
    )
end

