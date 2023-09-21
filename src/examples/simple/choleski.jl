using DecomposingPolynomialSystems
@var k[1:5]
K = [k[1] k[2] k[3]; 0 k[4] k[5]; 0 0 1]
F = System([sum(K'*K)]; variables=k)
A, S, _ = scalingSymmetries(F)

A