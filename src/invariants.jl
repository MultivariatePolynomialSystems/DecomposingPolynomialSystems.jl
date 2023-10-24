export invariants

function num_constraints(block_partition::Vector{Vector{Int}})::Int
    return length(block_partition)*(length(block_partition[1])-1)
end

function sumFromTo(m, n)
    return Int(n*(n+1)/2 - m*(m+1)/2)
end

function num_instances(n_unknowns::Int, n_params::Int, degree::Int, block_partition::Vector{Vector{Int}})::Int
    n_constraints = num_constraints(block_partition)
    n_all_mons = binomial(n_unknowns + n_params + degree, degree)
    n_param_mons = binomial(n_params + degree, degree)
    n_minors = sumFromTo(n_param_mons - 1, n_all_mons - 1)
    return Int(ceil(n_minors/n_constraints))
end

function create_vandermonde_matrix(block_partition::Vector{Vector{Int}}, eval_mons::Array{ComplexF64, 3}, eval_param_mons::Array{ComplexF64, 2})::Matrix{ComplexF64}
    n_instances = size(eval_mons, 3)
    n_blocks = length(block_partition)
    block_size = length(block_partition[1])
    n_constraints = num_constraints(block_partition)
    n_unknown_mons = size(eval_mons, 1)
    n_param_mons = size(eval_param_mons, 1)
    n_minors = sumFromTo(n_param_mons - 1, n_unknown_mons + n_param_mons - 1)
    
    A = zeros(ComplexF64, n_instances*n_constraints, n_minors)
    for i in 1:n_instances
        for j in 1:n_blocks
            block = block_partition[j]
            for k in 1:block_size-1
                row_idx = (i - 1)*n_constraints + (j-1)*(block_size-1) + k
                col_idx = 1
                for m in 1:n_unknown_mons
                    for n in (m+1):n_unknown_mons
                        M = eval_mons[[m,n], [block[k],block[k+1]], i]
                        A[row_idx, col_idx] = det(M)
                        col_idx += 1
                    end
                    for n in 1:n_param_mons
                        M = eval_mons[m, [block[k],block[k+1]], i]
                        A[row_idx, col_idx] = (M[1]-M[2])*eval_param_mons[n, i]
                        col_idx += 1
                    end
                end
            end
        end
    end

    return A
end

function interpolate_invariants(block_partition::Vector{Vector{Int}}, mons::Vector{Expression}, eval_mons::Array{ComplexF64, 3}, eval_param_mons::Array{ComplexF64, 2}, tol::Float64)::Matrix{ComplexF64}
    n_instances = size(eval_mons, 3)
    n_constraints = num_constraints(block_partition)
    n_unknown_mons = size(eval_mons, 1)
    n_param_mons = size(eval_param_mons, 1)
    n_minors = sumFromTo(n_param_mons - 1, n_unknown_mons + n_param_mons - 1)
    println("n_minors = ", n_minors)
    println("n_param_mons = ", n_param_mons)
    println("n_unknown_mons = ", n_unknown_mons)

    println("Creating vandermonde matrix of size ", (n_instances*n_constraints, n_minors), " from solutions...\n")
    A = create_vandermonde_matrix(block_partition, eval_mons, eval_param_mons)
    @assert size(A, 1) >= size(A, 2)

    println("Computing nullspace...")
    inv_minors = Matrix{ComplexF64}(transpose(nullspace(A)))
    println("Size of the transposed nullspace: ", size(inv_minors), "\n")

    if size(inv_minors, 1) != 0
        println("Computing reduced row echelon form of the transposed nullspace...")
        inv_minors = sparsify(rref(inv_minors, tol), tol, digits=1)
    end

    return inv_minors
end

function invariants(F::SampledSystem; degree::Int, tol::Float64=1e-5, param_dep::Bool=true)::Vector{Matrix{ComplexF64}}
    params = parameters(F.equations)
    param_dep ? vars = vcat(variables(F.equations), params) : vars = variables(F.equations)
    param_dep ? n_params = length(params) : n_params = 0
    
    mons = get_monomials_factorization(vars, degree, n_params)
    println(mons)
    println("Evaluating monomials...\n\n")
    eval_mons = evaluate_monomials_at_samples(mons, F.solutions, F.parameters, vars)
    println("eval_mons size = ", size(eval_mons))

    param_mons = get_monomials(params, degree)
    eval_param_mons = evaluate_monomials_at_samples(param_mons, F.parameters, params)

    invs_minors = Vector{Matrix{ComplexF64}}([])
    for i in 1:length(F.block_partitions)
        printstyled("Interpolating invariants for the ", i, "-th block partition...\n", color=:blue)
        push!(invs_minors, interpolate_invariants(F.block_partitions[i], mons, eval_mons, eval_param_mons, tol))
        println()
    end
    return invs_minors
end

function invariants(F::System, xp0::Tuple{Vector{ComplexF64}, Vector{ComplexF64}}; degree::Int=0, tol::Float64=1e-5, expected_n_sols::Int=0, param_dep::Bool=true)::Tuple{SampledSystem, Vector{FactorizingMap}}
    F = run_monodromy(F, xp0, expected_n_sols=expected_n_sols)

    println("\nNumber of nontrivial block partitions: ", length(F.block_partitions), "\n")
    if isempty(F.block_partitions)
        return (F, [])
    end

    n_unknowns = length(xp0[1])
    n_params = length(xp0[2])
    n_constraints = min([num_constraints(F.block_partitions[i]) for i in 1:length(F.block_partitions)]...)

    param_dep ? n_vars = n_unknowns + n_params : n_vars = n_unknowns
    n_instances = max(Int(ceil(1/n_constraints*binomial(n_vars + degree, degree))), binomial(n_params + degree - 1, degree - 1))

    F = sample_system(F, n_instances)

    return (F, compute_invariants(F, degree=degree, tol=tol, param_dep=param_dep))
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
    M = Matrix{Union{Nothing, Expression}}(hcat([[nothing for i in 1:n_sols] for j in 1:n_sols]...))
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
    return sparsify(vec(p2a(c)), 1e-5)
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