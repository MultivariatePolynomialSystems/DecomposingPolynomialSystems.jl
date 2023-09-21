using DecomposingPolynomialSystems
using LinearAlgebra

function fabricateSolution()
    R = c2R(randn(ComplexF64, 3, 1))
    t = randn(ComplexF64, 3)
    X = randn(ComplexF64, 3, 3)
    x = R*X + t*ones(1, 3)
    return (Vector{ComplexF64}(vcat(M2V(R), t, ones(3))), vcat(M2V(X), M2V(x)))
end

function p3p()
    @var R[1:3,1:3], t[1:3], α[1:3], X[1:3,1:3], x[1:3,1:3]
    eqs = vcat(reshape(transpose(R)*R - I, 9), [det(R) - 1])
    for i in 1:3
        append!(eqs, α[i]*x[:,i] - R*X[:,i] - t)
    end
    return System(eqs; variables = vcat(M2V(R), t, α), parameters = vcat(M2V(X), M2V(x)))
end


F = p3p()
xp0 = fabricateSolution()
norm(F(xp0[1], xp0[2]))


U, S, _ = scalingSymmetries(F)
U
S

symmetries = compute_symmetries(F, degree=1, param_dep=true)
symmetries[1]