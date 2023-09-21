using DecomposingPolynomialSystems, HomotopyContinuation

@var E[1:3,1:3], R[1:3,1:3], t[1:3]
eqs = vcat(M2V(E-xx(t)*R), M2V(R'*R-eye(3)), det(R)-1)
F = System(eqs; variables = vcat(M2V(E), M2V(R), t))

A, S = scalingSymmetries(F)

A = A[:,1:9]

MDs = multidegreesUpToTotalDegree(9, 3)
classes = partitionMultidegrees(MDs, A, S)
max(length.(classes)...)

v = M2V(2*E*E'*E - tr(E'*E)*E)
exponents_coefficients(v[4], M2V(E))[2]