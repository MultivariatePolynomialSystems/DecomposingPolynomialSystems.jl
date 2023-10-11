using DecomposingPolynomialSystems, HomotopyContinuation, LinearAlgebra

v = 2 # number of views
@var R[1:3,1:3,1:v], t[1:3,1:v], α[1:5,1:v], x[1:3,1:5,1:v], a[1:4]
eqs = vcat(reshape(transpose(R)*R - I, 9), [det(R) - 1])
for i in 1:5
    append!(eqs, β[i]*[y[:,i]; 1] - R*α[i]*[x[:,i]; 1] - t)
end
#push!(eqs, transpose(a)*[t; 1])

function fabricateSample()
    R₁ = eye(3)
    t₁ = zeros(CC, 3)
    R₂ = c2R(randn(CC, 3, 1))
    t₂ = randn(CC, 3)
    X = randn(CC, 3, 5)
    x = R₁*X + t₁*ones(1,5)
    y = R₂*X + t₂*ones(1,5)
    α = reshape(x[3,:], 1, 5)
    β = reshape(y[3,:], 1, 5)
    x = x./α
    y = y./β
    n = nullspace(reshape(vcat(t₂, 1), 1, 4))
    a = M2V(n*randn(CC, size(n, 2), 1))
    return (Vector{ComplexF64}(vcat(M2V(R₂), t₂, M2V(α), M2V(β))), vcat(M2V(x[1:2,:]), M2V(y[1:2,:]), a))
end

function eqs_5pp()
    
    return System(eqs; variables = vcat(reshape(R, 9), t, α, β), parameters = vcat(reshape(x, 10), reshape(y, 10)))
end

