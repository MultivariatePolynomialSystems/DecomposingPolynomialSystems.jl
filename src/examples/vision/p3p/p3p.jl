using DecomposingPolynomialSystems, LinearAlgebra

function fabricate_sample()
    R = c2R(randn(CC, 3, 1))
    t = randn(CC, 3)
    X = randn(CC, 3, 3)
    x = [R t]*a2p(X)
    α = ones(CC, 3)
    return ([M2V(R); t; α], [M2V(X); M2V(x)])
end

@var R[1:3,1:3], t[1:3], α[1:3], X[1:3,1:3], x[1:3,1:3]
eqs = vcat(reshape(transpose(R)*R - I, 9), [det(R) - 1])
for i in 1:3
    append!(eqs, α[i]*x[:,i] - [R t]*[X[:,i]; 1])
end
F = System(eqs; variables = vcat(M2V(R), t, α), parameters = vcat(M2V(X), M2V(x)))

F = run_monodromy(F)
scaling_symmetries(F)

symmetries_fixing_parameters!(F; degree_bound=1, param_dep=false)