using HomotopyContinuation, DecomposingPolynomialSystems, LinearAlgebra

@var fx fy s cx cy q[1:4,1:3] t[1:3, 1:3]
R1 = q2R(q[1:4,1])
t1 = t[:,1]
R2 = q2R(q[1:4,2])
t2 = t[:,2]
R3 = q2R(q[1:4,3])
t3 = t[:,3]
K = [fx s cx; 0 fy cy; 0 0 1]
Kinv = [fx 0 0; s fy 0; cx cy 1]
F12 = Kinv'*R2*xx(t2-t1)*R1'*Kinv
F13 = Kinv'*R3*xx(t3-t1)*R1'*Kinv
F23 = Kinv'*R3*xx(t3-t2)*R2'*Kinv

vars = vcat([fx, fy, s, cx, cy], M2V(q), M2V(t))

Jac_Fx = Matrix{Float64}(subs(differentiate(vcat(M2V(F12), M2V(F13), M2V(F23)), vars), vars => randn(length(vars))))
Ker_dFx = nullspace(Jac_Fx)

dim_im = length(vars) - size(Ker_dFx, 2)
