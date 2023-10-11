using DecomposingPolynomialSystems, AbstractAlgebra
@var R[1:3,1:3]
F = System(vcat(M2V(R'*R-eye(3)), [det(R)-1]))
A, S, K = scalingSymmetries(F, in_hnf=false)

hnf(matrix(GF(2), Matrix(K)))

A

S, T, U = snf_with_transform(K)
v*K

A = zeros(Int, 6, 9)
A[1,:] = [1, 1, 0, 1, 1, 0, 1, 1, 0]
A[2,:] = [0, 1, 1, 0, 1, 1, 0, 1, 1]
A[3,:] = [1, 0, 1, 1, 0, 1, 1, 0, 1]
A[4,:] = [1, 1, 1, 1, 1, 1, 0, 0, 0]
A[5,:] = [1, 1, 1, 0, 0, 0, 1, 1, 1]
A[6,:] = [0, 0, 0, 1, 1, 1, 1, 1, 1]

hnf(matrix(GF(2), A))