using HomotopyContinuation, DecomposingPolynomialSystems, LinearAlgebra

@var fx fy s cx cy H[1:4,1:4] q[1:4,1:3] t[1:3, 1:3]
R1 = q2R(q[1:4,1])
t1 = t[:,1]
R2 = q2R(q[1:4,2])
t2 = t[:,2]
R3 = q2R(q[1:4,3])
t3 = t[:,3]
K = [fx s cx; 0 fy cy; 0 0 1]
P1 = K*[R1 t1]*H
P2 = K*[R2 t2]*H
P3 = K*[R3 t3]*H

vars = vcat([fx, fy, s, cx, cy], M2V(H), M2V(q), M2V(t))

Jac_Fx = Matrix{Float64}(subs(differentiate(vcat(M2V(P1), M2V(P2), M2V(P3)), vars), vars => randn(length(vars))))
Ker_dFx = nullspace(Jac_Fx)

dim_im = length(vars) - size(Ker_dFx, 2)

using AbstractAlgebra
A = matrix(ZZ, 
    [1 1 0 0 1 1 0 0 -1 0 0 0 0; 
     0 0 1 1 0 0 1 1 0 -1 0 0 0; 
     1 1 1 1 0 0 0 0 0 0 -1 0 0; 
     0 0 0 0 1 1 1 1 0 0 0 -1 0;
     1 1 1 1 1 1 1 1 0 0 0 0 -1])
H = hnf(A)

