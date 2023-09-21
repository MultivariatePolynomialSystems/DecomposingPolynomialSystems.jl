# using HomotopyContinuation, AbstractAlgebra, Groebner, Combinatorics
# include("../../../../src/utils/utils_cpu.jl")
#
# # Getting the constraints by elimination... 24 hours isn't enough :-(
# @var P[1:2,1:n,1:n] T[1:2,1:2,1:2,1:2]
# str_vars = [string(var) for var in [vec(reshape(P, 1, 2*4*4)); vec(reshape(T, 1, 2^4))]]
#
# R, vars = PolynomialRing(QQ, str_vars, ordering=:lex)
#
# P, T = reshape(vars[1:2*4*4],2,4,4), reshape(vars[2*4*4+1:2*4*4+2^4],2,2,2,2)
#
# eqs = [T[i1,i2,i3,i4] - (-1)^(i1+i2+i3+i4)*exprDet([transpose(P[i1,:,1]); transpose(P[i2,:,2]); transpose(P[i3,:,3]); transpose(P[i4,:,4])], expnd=false) for (i1,i2,i3,i4) in vec(reshape(collect(Iterators.product(ntuple(_ -> 1:2, 4)...)), 2^4, 1))]
# eqs
#
# gb = groebner(eqs)
# using FileIO, JLD2
# FileIO.save("radial_quad_tensor_constraints_gb.jld2", "gb", gb)
#
#
#
# # Interpolating the constraints on CPU (nothing found up to degree 4)
# function Ps2T(Ps)
#     return [(-1)^(i1+i2+i3+i4)*det([transpose(Ps[1][i1,:]); transpose(Ps[2][i2,:]); transpose(Ps[3][i3,:]); transpose(Ps[4][i4,:])]) for (i1,i2,i3,i4) in vec(reshape(collect(Iterators.product(ntuple(_ -> 1:2, 4)...)), 2^4, 1))]
# end
#
# @var T[1:2,1:2,1:2,1:2]
# d = 4
# mons = get_monomials(reshape(T, 16), d)
# n_mons = length(mons)
#
# Ts = zeros(ComplexF64, n_mons+1, 16)
# for i=1:size(Ts,1)
#     Ps = [rand(ComplexF64, 2, 4) for i=1:4]
#     Ts[i,:] = Ps2T(Ps)
# end
#
# subs(mons, reshape(T, 16) => Ts[1,:])
#
# vandermondeMatrix = zeros(ComplexF64, size(Ts,1), n_mons)
# for i=1:size(vandermondeMatrix,1)
#     vandermondeMatrix[i,:] = Vector{ComplexF64}(subs(mons, reshape(T, 16) => Ts[i,:]))
# end
#
# nullspace(vandermondeMatrix)
#
#
#
# # Interpolating the constraints on GPU
# using CUDA
# using HomotopyContinuation
# include("utils/utils_cpu.jl")
# include("utils/utils_gpu.jl")
# function Ps2T(Ps)
#     return [(-1)^(i1+i2+i3+i4)*det([transpose(Ps[1][i1,:]); transpose(Ps[2][i2,:]); transpose(Ps[3][i3,:]); transpose(Ps[4][i4,:])]) for (i1,i2,i3,i4) in vec(reshape(collect(Iterators.product(ntuple(_ -> 1:2, 4)...)), 2^4, 1))]
# end
#
# @var T[1:2,1:2,1:2,1:2]
# d = 12
# mons = get_monomials_RQT_gpu(reshape(T, 16), d)
# n_mons = length(mons)
#
# Ts = CUDA.zeros(ComplexF64, n_mons+1, 16)
# for i=1:size(Ts,1)
#     Ps = [rand(ComplexF64, 2, 4) for i=1:4]
#     copyto!(view(Ts, i, :), CuArray(Ps2T(Ps)))
# end
#
# vandermondeMatrix = CUDA.ones(ComplexF64, size(Ts,1), n_mons)
# for i=1:size(vandermondeMatrix,1)
#     evaluate_monomials_RQT_gpu!(view(Ts,i,:), d, view(vandermondeMatrix, i, :))
# end
#
# tol = 1e-5
# nullspace_gpu(vandermondeMatrix, tol)
#
# ###############
# ###############
# ###############
using Oscar, HomotopyContinuation
include("../../../../src/utils/utils_cpu.jl")

@var P[1:4,1:4], w, T[1:16]
str_vars = [string(var) for var in vcat(reshape(P,16),[w],T)]

R, vars = PolynomialRing(GF(101), str_vars)
P, w, T = reshape(vars[1:16],4,4), vars[17], vars[18:33]

p11, p22, p33, p44 = -T[2], -T[3], -T[4], -T[5]
p21 = (-T[6]+p11*p22, p11)
p31 = (-T[7]+p11*p33, p11)
p41 = (-T[8]+p11*p44, p11)
p23p32 = -T[9]+p22*p33
p24p42 = -T[10]+p22*p44
p34p43 = -T[11]+p33*p44

eq1 = p11*T[12] - p11*p22*p31[1] + p11*P[2,3]*p31[1] + p11*p21[1]*P[3,2] - p11^2*p23p32 - p11*p21[1]*p33 + p11^2*p22*p33
eq2 = p11*T[13] - p11*p22*p41[1] + p11*P[2,4]*p41[1] + p11*p21[1]*P[4,2] - p11^2*p24p42 - p11*p21[1]*p44 + p11^2*p22*p44
eq3 = p11*T[14] - p11*p33*p41[1] + p11*P[3,4]*p41[1] + p11*p31[1]*P[4,3] - p11^2*p34p43 - p11*p31[1]*p44 + p11^2*p33*p44
eq4 = T[15] - p33*p24p42 + P[2,3]*P[3,4]*P[4,2] + P[2,4]*P[3,2]*P[4,3] - p22*p34p43 - p23p32*p44 + p22*p33*p44
A = p21[1]*T[11] - p31[1]*(P[2,3]*p44-P[4,3]*P[2,4]) + p41[1]*(P[2,3]*P[3,4]-p33*P[2,4])
B = p21[1]*(P[3,2]*p44-P[3,4]*P[4,2]) - p31[1]*T[10] + p41[1]*(p22*P[3,4]-P[3,2]*P[2,4])
C = p21[1]*(P[3,2]*P[4,3]-p33*P[4,2]) - p31[1]*(p22*P[4,3]-P[4,2]*P[2,3]) + p41[1]*T[9]
eq5 = T[16] - (p11*T[15] - A + B - C)
eq6 = -w*p11 + 1
eqs = [eq1, eq2, eq3, eq4, eq5, eq6]

I = ideal(eqs)

G = f4(I, eliminate=17)
