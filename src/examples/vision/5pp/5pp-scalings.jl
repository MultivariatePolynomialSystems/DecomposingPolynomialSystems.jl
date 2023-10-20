using DecomposingPolynomialSystems

@var R[1:3,1:3] t[1:3] α[1:5] β[1:5] x[1:3,1:5] y[1:3,1:5] a[1:4]
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
    a = vec(n*randn(CC, size(n, 2)))
    return (vcat(R₂[:], t₂, α, β), vcat(x[:], y[:], a))
end

xp0 = fabricateSample()

F = run_monodromy(F, xp0)
scalings = scaling_symmetries(F)
d = 3
MDs = multidegrees_up_to_total_degree(n_variables(F), d)
classes = partition_multidegrees(MDs, grading)

symmetries = symmetries_fixing_parameters_graded!(F, grading, MDs, classes)
@time symmetries = symmetries_fixing_parameters_graded!(F, grading, MDs, classes);
@profview symmetries = symmetries_fixing_parameters_graded!(F, grading, MDs, classes);

symmetries = symmetries_fixing_parameters!(F; degree=3)
