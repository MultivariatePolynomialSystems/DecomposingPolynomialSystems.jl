using GroebnerBasis, Singular

R, vars = PolynomialRing(QQ, ["w", "x", "a", "b"])
w,x,a,b = vars
I = Ideal(R, [a*x^4+b*x^2+a,w*x-1])
J = f4(I, monorder=:lex, reducegb=1)


using GroebnerBasis, Singular, HomotopyContinuation
include("../../utils/utils_cpu.jl")

@var Rot[1:3,1:3], t[1:3], E[1:3,1:3]
str_vars = [string(var) for var in vcat(reshape(Rot, 9), t, reshape(E, 9))]

R, vars = PolynomialRing(QQ, str_vars)
Rot, t, E = reshape(vars[1:9],3,3), vars[10:12], reshape(vars[13:end],3,3)

rotEqs = vcat(vcat(M2L(Matrix{}(transpose(Rot))*Rot-[1 0 0; 0 1 0; 0 0 1])...), [exprDet(Rot, expnd=false)-1])
essEqs = vcat(M2L(E - v_x(t)*Rot)...)

eqs = vcat(rotEqs, essEqs)
J = f4(Ideal(R, eqs), reducegb=1)
gens(J)


# using Oscar, HomotopyContinuation
# @var R[1:3,1:3], t[1:3], E[1:3,1:3]
# str_vars = [string(var) for var in vcat(reshape(R, 9), t, reshape(E, 9))]
#
# R, vars = PolynomialRing(GF(30011), str_vars)
#
# R, t, E = reshape(vars[1:9],3,3), vars[10:12], reshape(vars[13:end],3,3)
#
# rotEqs = vcat(vcat(M2L(Matrix{}(transpose(R))*R-[1 0 0; 0 1 0; 0 0 1])...), [exprDet(R, expnd=false)-1])
# essEqs = vcat(M2L(E - v_x(t)*R)...)
#
# eqs = vcat(rotEqs, essEqs)
#
# @time f4(ideal(eqs), eliminate=12)
