using HomotopyContinuation, GAP, LinearAlgebra
include("../../../../src/common.jl")
include("../../../../src/implicitization/implicitize.jl")
include("../../../../src/factorizing_maps/compute_factorizing_maps.jl")
include("../eqs_wo_depths.jl")

F = eqs_5pp()
xz0 = (x0, z0) = fabricateSolution()
norm(F(x0, z0))

F = run_monodromy(F, xz0)
F = sample_system(F, 100)

unknowns = variables(F.equations)
params = parameters(F.equations)
vars = vcat(unknowns, params)
mons = get_monomials_factorization(vars, 2, length(params))
evaluated_mons = evaluate_monomials_at_samples(mons, F.solutions, F.parameters, vars)
block_partition = F.block_partitions[1]

function create_vandermonde_matrix(block_partition::Vector{Vector{Int}}, evaluated_mons::Array{ComplexF64, 3})::Matrix{ComplexF64}
    n_constraints = length(block_partition)*(length(block_partition[1])-1)
    n_mons = size(evaluated_mons, 1)
    n_blocks = length(block_partition)
    block_size = length(block_partition[1])
    n_instances = size(evaluated_mons, 3)


    A = zeros(ComplexF64, n_instances*n_constraints, n_mons)
    @assert size(A, 1) >= size(A, 2)
    for i in 1:n_instances
        for j in 1:n_blocks
            block = block_partition[j]
            for k in 1:block_size-1
                row_idx = (i - 1)*n_constraints + (j-1)*(block_size-1) + k
                A[row_idx, :] = evaluated_mons[:, block[k], i] + evaluated_mons[:, block[k+1], i] # +/-
            end
        end
    end
    return A
end

A = create_vandermonde_matrix(block_partition, evaluated_mons)
φ_coeffs = Matrix{ComplexF64}(transpose(nullspace(A)))
tol = 1e-5
φ_coeffs = sparsify(rref(φ_coeffs, tol), tol, digits=1)

φ = Vector{Expression}([dot(φ_coeffs[i, :], mons) for i in 1:size(φ_coeffs, 1)])

for i = 1:length(φ)
    println(i, "-th function:")
    println(φ[i], "\n\n")
end

function is_param_dep_only(f::Expression, F::SampledSystem, tol::Float64)::Bool
    n_sols, n_instances = size(F.solutions, 2), size(F.solutions, 3)
    evals = zeros(ComplexF64, n_sols, n_instances)
    for i in 1:n_sols
        for j in 1:n_instances
            evals[i,j] = ComplexF64.(expand.(subs(f, vcat(F.equations.variables, F.equations.parameters) => vcat(F.solutions[:,i,j], F.parameters[:,j]))))
        end
    end
    return rank(evals, atol=tol) <= 1
end

function remove_param_dep_only(φ::Vector{Expression}, F::SampledSystem, tol::Float64)::Vector{Expression}
    nonparam_dep_only_φ = Vector{Expression}([])
    for i in 1:length(φ)
        if !is_param_dep_only(φ[i], F, tol)
            push!(nonparam_dep_only_φ, φ[i])
        end
    end
    return nonparam_dep_only_φ
end


tol = 1e-5
φ = remove_param_dep_only(φ, F, tol)
# Also removing φ[i] that are polynomials in the others φ[j] with parameter coefficients won't encrease template size.
φ = FactorizingMap(φ)

@var e[1:3,1:3]
e = M2V(e)
e_mons = get_monomials(e, 3)
exy_mons = vcat([[e[i]*params[j₁+2*k]*params[10+j₂+2*k] for i in 1:9 for j₁ in 1:2 for j₂ in 1:2] for k in 0:4]...)
ex_mons = [e[i]*params[j] for i in 1:9 for j in 1:10]
ey_mons = [e[i]*params[10+j] for i in 1:9 for j in 1:10]
new_mons = vcat(exy_mons, ex_mons, ey_mons, e_mons)

new_mons = get_monomials(vcat(e,params), 3)

G = implicitize(F, φ, block_partition, new_vars=e, mons=new_mons)
G.equations.expressions
for i = 1:length(G.equations.expressions)
    println(i, "-th function:")
    println(G.equations.expressions[i], "\n\n")
end
