using DecomposingPolynomialSystems

@var k[1:5] R[1:3,1:3,1:3] t[1:3,1:3] α[1:5,1:3] x[1:3,1:5,1:3] X[1:4,1:5]
K = [k[1] k[2] k[3]; 0 k[4] k[5]; 0 0 1]
Rs = [R[1:3,1:3,1], R[1:3,1:3,2], R[1:3,1:3,3]]
ts = [t[1:3,1], t[1:3,2], t[1:3,3]]
Ps = [K*[Rs[i] ts[i]] for i in 1:3]

rotEqs = vcat([vcat(M2V(R[:,:,i]'*R[:,:,i]-eye(3)), [det(R[:,:,i])-1]) for i in 1:3]...)
eqs = vcat(rotEqs, vcat([α[i,j]*x[:,i,j]-Ps[j]*X[:,i] for i in 1:5 for j in 1:3]...))
F = System(eqs)

A, S = scalingSymmetries(F)
S

A, S = scalingSymmetries(F, vcat(k, M2V(α), M2V(x)))
S
A