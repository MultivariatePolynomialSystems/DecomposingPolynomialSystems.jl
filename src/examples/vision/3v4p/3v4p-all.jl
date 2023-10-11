using DecomposingPolynomialSystems, HomotopyContinuation, LinearAlgebra

@var R[1:3,1:3,1:2] t[1:3,1:2] α[1:4,1:3] x[1:3,1:4,1:3]
eqs = vcat([vcat(M2V(R[:,:,i]'*R[:,:,i] - I), [det(R[:,:,i]) - 1]) for i in 1:2]...)
for i in 1:4
    for j in 1:2
        append!(eqs, α[i,j+1]*x[:,i,j+1] - R[:,:,j]*α[i,1]*x[:,i,1] - t[:,j])
    end
end
F = System(eqs; variables = vcat(M2V(R), M2V(t), M2V(α)), parameters = M2V(x))

A, S = scalingSymmetries(F)
