using DecomposingPolynomialSystems
using LinearAlgebra

function fabricateSolution()
    R1 = [1 0 0; 0 1 0; 0 0 1]
    t1 = [0; 0; 0]
    R2 = c2R(rand(ComplexF64, 3, 1))
    t2 = randn(ComplexF64, 3)
    X = randn(ComplexF64, 3, 5)
    x = R1*X + t1*ones(1,5)
    y = R2*X + t2*ones(1,5)
    x = x./reshape(x[3,:], 1, 5)
    y = y./reshape(y[3,:], 1, 5)
    n = nullspace(reshape(vcat(t2, 1), 1, 4))
    a = reshape(n*randn(ComplexF64, size(n, 2), 1), 4)
    return (Vector{ComplexF64}(vcat(reshape(R2, 9), t2)), vcat(reshape(x[1:2,:], 10), reshape(y[1:2,:], 10), a))
end

function eqs_5pp()
    @var R[1:3,1:3], t[1:3], x[1:2,1:5], y[1:2,1:5], a[1:4]
    eqs = vcat(reshape(transpose(R)*R - I, 9), [det(R) - 1])
    for i in 1:5
        push!(eqs, transpose([y[:,i]; 1])*xx(t)*R*[x[:,i]; 1])
    end
    push!(eqs, transpose(a)*[t; 1])
    return System(eqs; variables = vcat(reshape(R, 9), t), parameters = vcat(reshape(x, 10), reshape(y, 10), a))
end


F = eqs_5pp()
xp0 = fabricateSolution()
norm(F(xp0[1], xp0[2]))
F = run_monodromy(F, xp0)
deck_transformations  = compute_deck_transformations!(F, degree=3, param_dep=false)

for i = 1:12
    println(i, "-th rational function:")
    println(deck_transformations[2][i], "\n\n")
end
