using DecomposingPolynomialSystems
using HomotopyContinuation

@var x a b
F = System([(x)^2 + (a)*(x) + (b)]; variables=[x], parameters=[a,b])
F = run_monodromy(F)
sample_system!(F, 20)

p0 = F.parameters[:, 1]
sols0 = DecomposingPolynomialSystems.M2VV(F.solutions[:, :, 1])
res = solve(F.equations, sols0, start_parameters = p0, target_parameters = [[735,100]])

sols = DecomposingPolynomialSystems.VV2M(solutions(res[1][1]))
t = 100

log(abs(sols[1]/sols[2]))/log(t)
log(abs(sols[2]/sols[1]))/log(t)

sols = HomotopyContinuation.solutions(res[1][1])

deck_transformations = compute_deck_transformations!(F, degree=2, param_dep=true)



function f(a,b,c)
    return b^2 - 4*a*c
end

t = 100000
w = [-1.2, 0.4, 3.7]
log(abs(f(5*(t.^w)...)))/log(t)
log(abs(f(t.^(-w)...)))/log(t)

using DecomposingPolynomialSystems, HomotopyContinuation
import DecomposingPolynomialSystems: V2M, M2V, eye
@var R[1:3,1:3]

eqs = expand.(M2V((R - eye(3))'*(R - eye(3)) - eye(3)))