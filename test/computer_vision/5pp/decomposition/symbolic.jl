include("symbolic_tools.jl")

FF = GF(101)
R = buildvars("R", 1:3, 1:3)
t = buildvars("t", 1:3)
E = buildvars("E", 1:3, 1:3)
# E = buildvars("E", 1:3)

S, vars = PolynomialRing(FF, vcat(R,t,E), ordering=:degrevlex)

R = reshape(vars[1:9], 3, 3)
t = vars[10:12]
E = reshape(vars[13:21], 3, 3)
# E = vars[13:15]

rotEqs = vcat([exprDet(R, expnd=false)-1], reshape(transpose(R)*R-eye(3), 9))
essEqs = reshape(E - xx(t)*R, 9)
# essEqs = E - transpose(R)*t
eqs = vcat(rotEqs, essEqs)
for i in 1:5
    push!(eqs, transpose(rand(FF, 3))*xx(t)*R*rand(FF, 3))
end
# push!(eqs, transpose(rand(FF, 4))*[t; 1])

J = ideal(S, eqs)

K = eliminate(J, vars[1:12])
f4(K)
ltI = leading_ideal(K, ordering=lex(gens(S)))
leading_ideal(ltI, ordering=degrevlex(gens(S)))
mons = ltI.gens.O
exps = Vector{Vector{Int}}([])
for mon in mons
    push!(exps, collect(exponent_vector(mon, 1)))
end
exps




# Turning twisted pair into linear action
include("symbolic_tools.jl")

FF = GF(101)
t = buildvars("t", 1:3)
R = buildvars("R", 1:3, 1:3)
y = buildvars("y", 1:9)

S, vars = PolynomialRing(FF, vcat(t,R,y), ordering=:degrevlex)
t = vars[1:3]
R = reshape(vars[4:12], 3, 3)
y = vars[13:21]
p = M2V((2*t*transpose(t) - (transpose(t)*t).*eye(3))*R)
q = transpose(t)*t
eqs = [q*y[i] - p[i] for i in 1:9]
J = ideal(S, eqs)
f4(J)
