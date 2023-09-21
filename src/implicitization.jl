function evaluate_map_at_samples(factorizing_map::FactorizingMap, block_partition::Vector{Vector{Int}}, F::SampledSystem)::Array{ComplexF64, 3}
    n_unknowns = length(factorizing_map.map)
    n_blocks = length(block_partition)
    n_instances = size(F.parameters, 2)

    new_solutions = zeros(ComplexF64, n_unknowns, n_blocks, n_instances)

    for i in 1:n_instances
        params = F.parameters[:, i]
        for j in 1:n_blocks
            sol_idx = block_partition[j][1]
            sol = F.solutions[:, sol_idx, i]
            new_solutions[:, j, i] = ComplexF64.(expand.(subs(factorizing_map.map, vcat(variables(F.equations), parameters(F.equations)) => vcat(sol, params))))
        end
    end

    return new_solutions
end

function implicitize(F::SampledSystem, factorizing_map::FactorizingMap, block_partition::Vector{Vector{Int}}; new_vars::Vector{Variable}=Vector{Variable}([]), mons::Vector{Expression}=Vector{Expression}([]), degree::Int64=0, tol::Float64=1e-5)::SampledSystem
    if isempty(factorizing_map.map)
        return SampledSystem()
    end

    new_solutions = evaluate_map_at_samples(factorizing_map, block_partition, F)

    if mons == []
        mons = get_monomials(vcat(new_vars, parameters(F.equations)), degree)
    end
    if new_vars == []
        @var y[1:length(factorizing_map.map)]
        new_vars = y
    end
    evaluated_mons = evaluate_monomials_at_samples(mons, new_solutions, F.parameters, vcat(new_vars, parameters(F.equations)))

    m, n, k = size(evaluated_mons)
    A = Matrix{ComplexF64}(transpose(reshape(evaluated_mons, m, n*k)))
    # A = A[1:m, 1:m]
    @assert size(A, 1) >= size(A, 2)

    coeffs = Matrix{ComplexF64}(transpose(nullspace(A)))

    if size(coeffs, 1) != 0
        println("Computing the reduced row echelon form of the transposed nullspace...")
        coeffs = sparsify(rref(coeffs, tol), tol, digits=1)
    end

    G = System([dot(coeffs[i,:], mons) for i in 1:size(coeffs, 1)], variables=new_vars, parameters=parameters(F.equations))

    mon = action_on_blocks(F.monodromy_group, block_partition)
    block_partitions = all_block_partitions(mon)
    deck = centralizer(mon)

    return SampledSystem(G, new_solutions, F.parameters, mon, block_partitions, deck)
end
