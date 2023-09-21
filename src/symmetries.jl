include("interpolation.jl")

export FiniteSymmetryGroup
export scaling_symmetries, scaling_symmetries!
export symmetries_fixing_parameters!, symmetries_fixing_parameters_graded!
export symmetries_fixing_parameters

struct FiniteSymmetryGroup
    maps::Vector{Vector{Union{Nothing, Expression}}}
    vars::Vector{Variable}
    structure::GapObj
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
    if s == []
        return [], []
    end
    idxs = [findall(x->x==el, S) for el in s]
    Us = []
    for i in eachindex(idxs)
        push!(Us, U[idxs[i], :])
    end
    return s, Us # TODO: what if max(Us[i]) > MAX_INT64?
end

function _hnf_reduce(grading::Grading)::Grading
    for (i, (sᵢ, Uᵢ)) in enumerate(grading)
        if sᵢ == 0
            grading[i] = (sᵢ, Matrix(hnf(matrix(ZZ, Uᵢ))))
        else
            grading[i] = (sᵢ, lift.(Matrix(hnf(matrix(GF(sᵢ), Uᵢ)))))
        end
    end
    return grading
end

function scaling_symmetries(F::System; in_hnf::Bool=true)::Grading
    s, U = _snf_scaling_symmetries(F)
    if s == []
        return []
    end
    grading = collect(zip(s, U))
    return in_hnf ? _hnf_reduce(grading) : grading
end

function scaling_symmetries!(F::SampledSystem; in_hnf::Bool=true)::Grading
    gr = scaling_symmetries(F.system; in_hnf=in_hnf)
    F.grading = copy(gr)
    return gr
end

function _remove_zero_rows(M::Matrix)::Matrix
    return transpose(VV2M(filter(!iszero, M2VV(transpose(M)))))
end

function _remove_zero_rows(grading::Grading)::Grading
    for (i, (sᵢ, Uᵢ)) in enumerate(grading)
        grading[i] = (sᵢ, _remove_zero_rows(Uᵢ))
    end
    return grading
end

function _remove_dependencies(grading::Grading)::Grading
    # TODO: remove rows dependent on other blocks
    return grading
end

# TODO: what if vars is not a subset of Fvars?
function scaling_symmetries(F::System, vars::Vector{Variable})::Grading
    Fvars = vcat(F.variables, F.parameters)
    idx = [findfirst(v->v==var, Fvars) for var in vars]
    s, U = _snf_scaling_symmetries(F)
    U = [u[:, idx] for u in U]
    grading = _remove_zero_rows(_hnf_reduce(collect(zip(s, U))))
    return _remove_dependencies(grading)
end

function scaling_symmetries!(F::SampledSystem, vars::Vector{Variable})::Grading
    gr = scaling_symmetries(F.system, vars)
    F.grading = copy(gr)
    return gr
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
function _remove_zero_nums_and_denoms(
    coeffs::Matrix{CC},
    num_mons::Vector{Expression},
    denom_mons::Vector{Expression};
    logging::Bool=false
)::Matrix{CC}

    reasonable_rows = []
    n_num_mons, n_denom_mons = length(num_mons), length(denom_mons)
    @assert size(coeffs, 2) == n_num_mons + n_denom_mons
    for i in 1:size(coeffs, 1)
        if (!all(iszero, coeffs[i, 1:n_num_mons]) && !all(iszero, coeffs[i, n_num_mons+1:end]))
            push!(reasonable_rows, i)
        elseif logging
            println("Removed: ",
                dot(coeffs[i, 1:n_num_mons], num_mons) + dot(coeffs[i, n_num_mons+1:end],
                denom_mons)
            )
        end
    end
    return coeffs[reasonable_rows, :]
end

function _remove_zero_nums_and_denoms(
    coeffs::Matrix{CC},
    mons::Vector{Expression};
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

function _all_interpolated(symmetry_group::FiniteSymmetryGroup)::Bool
    all_interpolated = true
    for symmetry in symmetry_group.maps
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

function _init_symmetries(n_symmetries::Int, unknowns::Vector{Variable})::FiniteSymmetryGroup
    symmetries = [[nothing for j in eachindex(unknowns)] for i in 1:n_symmetries]
    symmetries = Vector{SymmetryMap}(symmetries)
    symmetries[1] = Expression.(unknowns)  # set the first to the identity
    return FiniteSymmetryGroup(symmetries, unknowns, GapObj())
end

function _print_symmetry_function(symmetries, symmetry_id)

end

function _interpolate_symmetry_function(
    permutation::Vector{Int},
    values::SubArray{CC, 2},
    eval_num_mons::Array{CC, 3},
    eval_denom_mons::Array{CC, 3},
    num_mons::Vector{Expression},
    denom_mons::Vector{Expression},
    tol::Float64;
    logging::Bool=false
)::Union{Nothing, Expression}

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
    return rational_function(coeffs, num_mons, denom_mons; printing=false)
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
    grading::Grading, 
    mds::Vector{Multidegree}, 
    classes::Dict{Vector{Int}, Vector{Int}}; 
    tol::Float64=1e-5,
)::Vector{SymmetryMap}

    unkns, vars = unknowns(F), variables(F)
    
    max_n_mons = max(length.(collect(values(classes)))...) # size of the largest class
    n_unknowns, n_sols, _ = size(F.samples.solutions)
    n_instances = Int(ceil(2/n_sols*max_n_mons))

    C = F.symmetry_permutations
    symmetries = _init_symmetries(length(C), unkns)

    sample_system!(F, n_instances)
    
    for (num_deg, num_ids) in classes
        num_mds = mds[num_ids]
        num_mons = mds2mons(num_mds, vars)
        eval_num_mons = nothing
        for i in 1:n_unknowns
            denom_deg = _num_deg2denom_deg(num_deg, grading, i) # i-th variable
            denom_ids = get(classes, denom_deg, nothing)
            if !isnothing(denom_ids)
                denom_mds = mds[denom_ids]
                g = gcd_mds([gcd_mds(num_mds), gcd_mds(denom_mds)])
                if iszero(g) & !only_param_dep(vcat(num_mds, denom_mds), n_unknowns)
                    if isnothing(eval_num_mons)
                        eval_num_mons = evaluate_monomials_at_samples_(num_mds, F.samples)
                    end
                    denom_mons = mds2mons(denom_mds, vars)
                    eval_denom_mons = evaluate_monomials_at_samples_(denom_mds, F.samples)
                    for j in 2:length(C)
                        symmetry = symmetries[j]
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
                                    unkns[i],
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
            return symmetries
        end
    end

    return symmetries
end

function symmetries_fixing_parameters_graded!(
    F::SampledSystem,
    grading::Grading;
    degree_bound::Int,
    tol::Float64=1e-5,
    param_dep::Bool=true
)::FiniteSymmetryGroup

    MDs = multidegrees_up_to_total_degree(n_variables(F), degree_bound)
    classes = partition_multidegrees(MDs, grading)
    return symmetries_fixing_parameters_graded!(
        F,
        grading,
        MDs,
        classes;
        tol=tol
    )
end

function symmetries_fixing_parameters!(
    F::SampledSystem; 
    degree_bound::Int, 
    tol::Float64=1e-5, 
    param_dep::Bool=true
)::Vector{SymmetryMap}

    n_unknowns, n_sols, _ = size(F.samples.solutions)
    param_dep ? vars = variables(F) : vars = unknowns(F)
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
            return symmetries
        end
    end

    return symmetries
end

function symmetries_fixing_parameters(
    F::System, 
    xp0::Tuple{Vector{CC}, Vector{CC}}; 
    degree_bound::Int, 
    tol::Float64=1e-5, 
    param_dep::Bool=true,
    monodromy_options::Tuple=()
)::Tuple{FiniteSymmetryGroup, SampledSystem}

    F = run_monodromy(F, xp0; monodromy_options...)

    println("Order of the symmetry group: ", length(F.symmetry_permutations), "\n")
    if length(F.symmetry_permutations) == 1 # trivial group of symmetries
        return [Expression.(unknowns(F))] # return the identity map
    end

    grading = scaling_symmetries!(F)
    if length(grading) == 0
        return (symmetries_fixing_parameters!(
            F;
            degree_bound=degree_bound,
            tol=tol,
            param_dep=param_dep
            ), F)
    else
        println("Found non-trivial grading:")
        for (sᵢ, Uᵢ) in grading
            printstyled(size(Uᵢ, 1), " scalings of order ", sᵢ, "\n", color=:green)
        end
        printstyled("Running graded version...\n", color=:green)
        return (symmetries_fixing_parameters_graded!(
            F,
            grading;
            degree_bound=degree_bound,
            tol=tol,
            param_dep=param_dep
            ), F)
    end
end

"""
    symmetries_fixing_parameters(F; degree)
This is an example of Docstring. This function returns the
group of symmetries of the polynomial system `F` that fix the 
parameters.
```math
Deck(f)
```
"""
function symmetries_fixing_parameters(
    F::System;
    degree_bound::Int,
    tol::Float64=1e-5,
    param_dep::Bool=true,
    monodromy_options::Tuple=()
)::Tuple{FiniteSymmetryGroup, SampledSystem}

    xp0 = HomotopyContinuation.find_start_pair(F)
    return symmetries_fixing_parameters(
        F,
        xp0;
        degree_bound=degree_bound,
        tol=tol,
        param_dep=param_dep,
        monodromy_options=monodromy_options
    )
end

