using DecomposingPolynomialSystems
@var x[1:2] p[1:2]
F = System([x[1]^2 - x[2]^2 - p[1], 2*x[1]*x[2] - p[2]]; variables=x, parameters=p)
F = run_monodromy(F)
deck = symmetries_fixing_parameters!(F; degree_bound=1, param_dep=false)
deck[1]

#-------------------------------------------------------------------------------------------
using DecomposingPolynomialSystems, HomotopyContinuation
@var x a b
F = System([x^2+a*x+b]; variables=[x], parameters=[a,b])
F = run_monodromy(F)
deck = symmetries_fixing_parameters!(F; degree_bound=1, param_dep=true)

scaling_symmetries(F)

#-------------------------------------------------------------------------------------------
using DecomposingPolynomialSystems
@var x y z w a b c
eqs = [2*w^2+1, 2*x+4*w*y+2*a, -3*y^2-z^2+4*a*w*y+2*b, -w*y^3-3*w*y*z^2-a*y^2-a*z^2+2*b*w*y+2*c]
F = System(eqs; variables=[x,y,z,w], parameters=[a,b,c])
F = run_monodromy(F)
deck = symmetries_fixing_parameters!(F; degree_bound=1, param_dep=false)

@var x a b
md = Int8[1,2,3]
vars = [x,a,b]
m = Monomial(md, vars)
