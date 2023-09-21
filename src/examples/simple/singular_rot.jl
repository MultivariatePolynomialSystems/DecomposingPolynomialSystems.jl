using DecomposingPolynomialSystems
using HomotopyContinuation
using LinearAlgebra

CC = ComplexF64

@var R[1:3,1:3]
F = System(M2V(R'*R-eye(3)))

R₀ = c2R(randn(CC, 3))
R₀ = eye(3)
svdvals(jac(F, CC.(M2V(R₀))))