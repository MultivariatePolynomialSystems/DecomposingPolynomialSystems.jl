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

d = Int8(3)
mons = monomials(variables(F), d)

evaluate(mons, F.samples)
@btime evaluate(mons, F.samples);

evaluate_monomials_at_samples_(MDs, F.samples)
@btime evaluate_monomials_at_samples_(MDs, F.samples);

vars = DecomposingPolynomialSystems.variables(F)
mons = mds2mons(MDs, vars)
evaluate_monomials_at_samples(mons, F.samples, vars) # 40 sec
@btime evaluate_monomials_at_samples(mons, F.samples, vars);
@profview evaluate_monomials_at_samples(mons, F.samples, vars);
