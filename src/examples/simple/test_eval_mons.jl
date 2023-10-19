using DecomposingPolynomialSystems
using BenchmarkTools

@var R[1:3,1:3] t[1:3] α[1:5], β[1:5], x[1:3,1:5], y[1:3,1:5], a[1:4]
eqs = vcat((R'*R - I)[:], [det(R) - 1])
for i in 1:5
    append!(eqs, β[i]*y[:,i] - R*α[i]*x[:,i] - t)
end
push!(eqs, a'*[t; 1])

F = System(eqs; variables = vcat(R[:], t, α, β), parameters = vcat(x[:], y[:], a))

function fabricateSample()
    R₁ = eye(CC, 3)
    t₁ = zeros(CC, 3)
    R₂ = c2R(randn(CC, 3))
    t₂ = randn(CC, 3)
    X = randn(CC, 3, 5)
    x = [R₁ t₁]*a2p(X)
    y = [R₂ t₂]*a2p(X)
    α, β = ones(CC, 5), ones(CC, 5)
    n = nullspace(reshape([t₂; 1], 1, 4))
    a = n*randn(CC, size(n, 2))[:]
    return (vcat(R₂[:], t₂, α, β), vcat(x[:], y[:], a))
end

xp0 = fabricateSample()
F = run_monodromy(F, xp0)
sample_system!(F, 10)

function evaluate1(mons::MonomialVector, samples::VarietySamples)
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
            evaluated_mons[j, :, i] = vec(prod(sols.^mds[j][1:n_unknowns], dims=1)).*params_eval[j]
        end
    end
    return evaluated_mons
end

function evaluate2(mons::MonomialVector, samples::VarietySamples)
    solutions = samples.solutions
    parameters = samples.parameters

    n_unknowns, n_sols, n_instances = size(solutions)
    mds = mons.mds
    n_mds = length(mds)

    evaluated_mons = zeros(CC, n_mds, n_sols, n_instances)
    for i in 1:n_instances
        params = view(parameters, :, i)
        params_eval = [prod(params.^md[n_unknowns+1:end]) for md in mds]
        sols = view(solutions, :, :, i)
        for j in 1:n_mds
            evaluated_mons[j, :, i] = vec(prod(sols.^mds[j][1:n_unknowns], dims=1)).*params_eval[j]
        end
    end
    return evaluated_mons
end

function evaluate3(mons::MonomialVector, samples::VarietySamples)
    solutions = samples.solutions
    parameters = samples.parameters

    n_unknowns, n_sols, n_instances = size(solutions)
    mds = mons.mds
    n_mds = length(mds)

    nonzero_ids_unknowns = [findall(!iszero, md[1:n_unknowns]) for md in mds]
    nonzero_ids_params = [findall(!iszero, md[n_unknowns+1:end]) for md in mds]

    evaluated_mons = zeros(CC, n_mds, n_sols, n_instances)
    for i in 1:n_instances
        params = view(parameters, :, i)
        sols = view(solutions, :, :, i)
        for (j, md) in enumerate(mds)
            params_part = params[nonzero_ids_params[j]]
            md_params_part = md[n_unknowns+1:end][nonzero_ids_params[j]]
            params_eval = prod(params_part.^md_params_part)
            sols_part = view(sols, nonzero_ids_unknowns[j], :)
            md_sols_part = md[1:n_unknowns][nonzero_ids_unknowns[j]]
            evaluated_mons[j, :, i] = vec(prod(sols_part.^md_sols_part, dims=1)).*params_eval
        end
    end
    return evaluated_mons
end

d = 2
monsInt8 = monomials(variables(F), Int8(d))
monsInt = monomials(variables(F), d)

@btime evaluate1(monsInt8, F.samples);
@btime evaluate1(monsInt, F.samples);

@btime evaluate2(monsInt8, F.samples);
@btime evaluate2(monsInt, F.samples);

# THE FASTEST
@btime evaluate3(monsInt8, F.samples);
@btime evaluate3(monsInt, F.samples);