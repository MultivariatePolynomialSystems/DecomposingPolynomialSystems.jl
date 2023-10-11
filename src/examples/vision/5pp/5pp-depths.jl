using DecomposingPolynomialSystems, HomotopyContinuation, LinearAlgebra

function fabricateSolution()
    R1 = eye(3)
    t1 = zeros(3)
    R2 = c2R(randn(ComplexF64, 3))
    t2 = randn(ComplexF64, 3)
    X = randn(ComplexF64, 3, 5)
    x = [R1 t1]*a2p(X)
    y = [R2 t2]*a2p(X)
    α = ones(4)
    β = ones(5)
    return (Vector{ComplexF64}(vcat(α, β)), vcat(M2V(x), M2V(y)))
end

function eqs_5pp()
    @var α[1:5], β[1:5], x[1:3,1:5], y[1:3,1:5]
    eqs = []
    for i in 1:5
        for j in i:5
            if i != j
                dA = α[i]*x[:,i]-α[j]*x[:,j]
                dB = β[i]*y[:,i]-β[j]*y[:,j]
                append!(eqs, transpose(dA)*dA - transpose(dB)*dB)
            end
        end
    end
    return System(eqs; variables = vcat(α, β), parameters = vcat(M2V(x), M2V(y)))
end

F = eqs_5pp()
xp0 = fabricateSolution()
norm(F(xp0[1], xp0[2]))
run_monodromy(F, xp0)

A, S = scalingSymmetries(F)
MDs = multidegreesUpToTotalDegree(40, 6)
classes = partitionMultidegrees(MDs, A, S)
max(length.(classes)...)

sample_system!(F, 200)
compute_invariants(F, degree=2)


B = F.block_partitions[1]
sols = F.solutions[:,:,1]

A = zeros(ComplexF64, 10, 2)
for i in 1:10
    A[i,:] = [sols[1,B[i][1]]^2 - sols[1,B[i][2]]^2, sols[1,B[i][1]] - sols[1,B[i][2]]]
end
nullspace(A)

vars = variables(F.equations)
α = [1; vars[1:4]]
β = vars[5:9]
e = [β[i]^2/α[i]^2 for i in 1:5]
eval_at_sols(F, e)
are_LI(F, e)

e_inv = [α[i]^2/β[i]^2 for i in 1:5]
eval_at_sols(F, e)

e = [β[1]^(2*i) for i in 0:9]
are_LI(F, e)

sample_system!(F, 300)

symmetries, Fs = compute_symmetries(F, xp0; degree=3, std_form=true)
symmetries.basis
symmetries.coefficients[1]
symmetries.coefficients[2]