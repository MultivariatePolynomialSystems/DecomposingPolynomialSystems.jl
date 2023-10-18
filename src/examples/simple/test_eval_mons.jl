using DecomposingPolynomialSystems
using HomotopyContinuation, LinearAlgebra, BenchmarkTools

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

d = Int8(3)
n_variables(F)
MDs = multidegrees_affine(n_variables(F), d)

evaluate_monomials_at_samples(MDs, F.samples) # < 2 sec
@btime evaluate_monomials_at_samples(MDs, F.samples);

evaluate_monomials_at_samples_(MDs, F.samples)
@btime evaluate_monomials_at_samples_(MDs, F.samples);

vars = DecomposingPolynomialSystems.variables(F)
mons = mds2mons(MDs, vars)
evaluate_monomials_at_samples(mons, F.samples, vars) # 40 sec
@btime evaluate_monomials_at_samples(mons, F.samples, vars);
@profview evaluate_monomials_at_samples(mons, F.samples, vars);
