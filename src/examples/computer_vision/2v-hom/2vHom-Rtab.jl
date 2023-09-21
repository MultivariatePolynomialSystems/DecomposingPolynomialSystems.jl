using DecomposingPolynomialSystems, HomotopyContinuation, LinearAlgebra

@var R[1:3,1:3] t[1:3] α[1:4] β[1:4] x[1:3,1:4] y[1:3,1:4]
eqs = vcat(reshape(transpose(R)*R - I, 9), [det(R) - 1])
for i in 1:4
    append!(eqs, β[i]*y[:,i] - R*α[i]*x[:,i] - t)
end
push!(eqs, det(VV2M([a2p(α[i]*x[:,i]) for i in 1:4])))
F = System(eqs; variables = vcat(M2V(R), t, α, β), parameters = vcat(M2V(x), M2V(y)))

A, S = scalingSymmetries(F)
A
S

MDs = multidegreesUpToTotalDegree(44, 9)
classes = partitionMultidegrees(MDs, A, S)
max(length.(classes)...)