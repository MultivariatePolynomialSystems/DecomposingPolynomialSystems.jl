include("../../../src/utils/utils_cpu.jl")
include("../../../src/deck_transformations/compute_deck_transformations.jl")

function fabricateSolution()
    R1 = [1 0 0; 0 1 0; 0 0 1]
    t1 = [0; 0; 0]
    R2 = cayley_rot(randn(ComplexF64, 3, 1))
    t2 = randn(ComplexF64, 3)
    X = randn(ComplexF64, 3, 5)
    x = R1*X + t1*ones(1,5)
    y = R2*X + t2*ones(1,5)
    α = reshape(x[3,:], 1, 5)
    β = reshape(y[3,:], 1, 5)
    x = x./α
    y = y./β
    n = nullspace(reshape(vcat(t2, 1), 1, 4))
    a = reshape(n*randn(ComplexF64, size(n, 2), 1), 4)
    return (Vector{ComplexF64}(vcat(reshape(R2, 9), t2, reshape(α, 5), reshape(β, 5))), vcat(reshape(x[1:2,:], 10), reshape(y[1:2,:], 10), a))
end

function eqs_5pp()
    @var R[1:3,1:3], t[1:3], α[1:5], β[1:5], x[1:2,1:5], y[1:2,1:5], a[1:4]
    eqs = vcat(reshape(transpose(R)*R - I, 9), [det(R) - 1])
    for i in 1:5
        append!(eqs, β[i]*[y[:,i]; 1] - R*α[i]*[x[:,i]; 1] - t)
    end
    push!(eqs, transpose(a)*[t; 1])
    return System(eqs; variables = vcat(reshape(R, 9), t, α, β), parameters = vcat(reshape(x, 10), reshape(y, 10), a))
end


F = eqs_5pp()
(x0, p0) = fabricateSolution()
norm(F(x0, p0))
