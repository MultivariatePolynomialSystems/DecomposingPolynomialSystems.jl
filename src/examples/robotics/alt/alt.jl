using DecomposingPolynomialSystems, HomotopyContinuation

@var x y a b x̂ ŷ â b̂
@var p[1:8], p̂[1:8]

f1 = (x̂-ŷ)*(x-y);
f2 = (x-y)*(â*x̂-2*â*ŷ+2*b̂*x̂-b̂*ŷ);
f3 = (â^2*ŷ-2*â*b̂*x̂+2*â*b̂*ŷ-b̂^2*x̂)*(x-y);
f4 = â*b̂*(x-y)*(â*ŷ-b̂*x̂);
f5 = (x̂-ŷ)*(a*x-2*a*y+2*b*x-b*y);
f6 = -a*â*x*x̂+a*â*x*ŷ+a*â*x̂*y-2*a*â*y*ŷ-2*a*b̂*x*x̂+a*b̂*x*ŷ+4*a*b̂*x̂*y-2*a*b̂*y*ŷ+a*x*x̂*ŷ+a*x̂^2*y-2*a*x̂*y*ŷ-2*â*b*x*x̂+4*â*b*x*ŷ+â*b*x̂*y-2*â*b*y*ŷ+â*x^2*ŷ+â*x*x̂*y-2*â*x*y*ŷ-2*b*b̂*x*x̂+b*b̂*x*ŷ+b*b̂*x̂*y-b*b̂*y*ŷ-2*b*x*x̂*ŷ+b*x*ŷ^2+b*x̂*y*ŷ-2*b̂*x*x̂*y+b̂*x*y*ŷ+b̂*x̂*y^2-x^2*ŷ^2+2*x*x̂*y*ŷ-x̂^2*y^2;
f7 = 2*a*â*b̂*x*x̂-a*â*b̂*x*ŷ-2*a*â*b̂*x̂*y+2*a*â*b̂*y*ŷ-a*â*x*x̂*ŷ+2*a*â*x̂*y*ŷ+a*b̂^2*x*x̂-2*a*b̂^2*x̂*y-a*b̂*x*x̂*ŷ-2*a*b̂*x̂^2*y+2*a*b̂*x̂*y*ŷ-2*â^2*b*x*ŷ+â^2*b*y*ŷ-â^2*x^2*ŷ+2*â^2*x*y*ŷ+2*â*b*b̂*x*x̂-2*â*b*b̂*x*ŷ-â*b*b̂*x̂*y+2*â*b*b̂*y*ŷ+2*â*b*x*x̂*ŷ-2*â*b*x*ŷ^2-â*b*x̂*y*ŷ-â*b̂*x^2*ŷ-â*b̂*x̂*y^2+2*â*x^2*ŷ^2-2*â*x*x̂*y*ŷ+2*b*b̂*x*x̂*ŷ-b*b̂*x̂*y*ŷ+2*b̂^2*x*x̂*y-b̂^2*x̂*y^2-2*b̂*x*x̂*y*ŷ+2*b̂*x̂^2*y^2;
f8 = -a*â*b̂^2*x*x̂+a*â*b̂^2*x̂*y+a*â*b̂*x*x̂*ŷ-2*a*â*b̂*x̂*y*ŷ+a*b̂^2*x̂^2*y+â^2*b*b̂*x*ŷ-â^2*b*b̂*y*ŷ+â^2*b*x*ŷ^2+â^2*b̂*x^2*ŷ-â^2*b̂*x*y*ŷ-â^2*x^2*ŷ^2-2*â*b*b̂*x*x̂*ŷ+â*b*b̂*x̂*y*ŷ-â*b̂^2*x*x̂*y+â*b̂^2*x̂*y^2+2*â*b̂*x*x̂*y*ŷ-b̂^2*x̂^2*y^2;
f9 = (x̂-ŷ)*(a^2*y-2*a*b*x+2*a*b*y-b^2*x);
f10 = -2*a^2*b̂*x̂*y+a^2*b̂*y*ŷ-a^2*x̂^2*y+2*a^2*x̂*y*ŷ+2*a*â*b*x*x̂-2*a*â*b*x*ŷ-a*â*b*x̂*y+2*a*â*b*y*ŷ-a*â*x*x̂*y+2*a*â*x*y*ŷ+2*a*b*b̂*x*x̂-a*b*b̂*x*ŷ-2*a*b*b̂*x̂*y+2*a*b*b̂*y*ŷ-a*b*x*ŷ^2-a*b*x̂^2*y+2*a*b̂*x*x̂*y-a*b̂*x*y*ŷ-2*a*b̂*x̂*y^2-2*a*x*x̂*y*ŷ+2*a*x̂^2*y^2+â*b^2*x*x̂-2*â*b^2*x*ŷ-2*â*b*x^2*ŷ-â*b*x*x̂*y+2*â*b*x*y*ŷ+2*b^2*x*x̂*ŷ-b^2*x*ŷ^2+2*b*b̂*x*x̂*y-b*b̂*x*y*ŷ+2*b*x^2*ŷ^2-2*b*x*x̂*y*ŷ;
f11 = a^2*b̂^2*x̂*y+2*a^2*b̂*x̂^2*y-2*a^2*b̂*x̂*y*ŷ-a^2*x̂^2*y*ŷ-2*a*â*b*b̂*x*x̂+a*â*b*b̂*x*ŷ+a*â*b*b̂*x̂*y-2*a*â*b*b̂*y*ŷ+a*â*b*x*ŷ^2-a*â*b*x̂*y*ŷ-a*â*b̂*x*y*ŷ+a*â*b̂*x̂*y^2-a*b*b̂*x*x̂*ŷ+a*b*b̂*x̂^2*y+a*b*x*x̂*ŷ^2+a*b*x̂^2*y*ŷ-2*a*b̂^2*x*x̂*y+2*a*b̂^2*x̂*y^2+3*a*b̂*x*x̂*y*ŷ-3*a*b̂*x̂^2*y^2+â^2*b^2*x*ŷ+2*â^2*b*x^2*ŷ-2*â^2*b*x*y*ŷ-â^2*x^2*y*ŷ-2*â*b^2*x*x̂*ŷ+2*â*b^2*x*ŷ^2+â*b*b̂*x^2*ŷ-â*b*b̂*x*x̂*y-3*â*b*x^2*ŷ^2+3*â*b*x*x̂*y*ŷ+â*b̂*x^2*y*ŷ+â*b̂*x*x̂*y^2-b^2*x*x̂*ŷ^2-b̂^2*x*x̂*y^2;
f12 = (a*b̂*x̂*y-â*b*x*ŷ)*(a*b̂*x̂-a*x̂*ŷ-â*b*ŷ-â*b̂*x+â*b̂*y+â*x*ŷ+b*x̂*ŷ-b̂*x̂*y);
f13 = a*b*(x̂-ŷ)*(a*y-b*x);
f14 = a^2*b*b̂*x̂*y-a^2*b*b̂*y*ŷ+a^2*b*x̂^2*y-a^2*b*x̂*y*ŷ+a^2*b̂*x̂*y^2-a^2*x̂^2*y^2-a*â*b^2*x*x̂+a*â*b^2*x*ŷ+a*â*b*x*x̂*y-2*a*â*b*x*y*ŷ-a*b^2*x*x̂*ŷ+a*b^2*x*ŷ^2-2*a*b*b̂*x*x̂*y+a*b*b̂*x*y*ŷ+2*a*b*x*x̂*y*ŷ+â*b^2*x^2*ŷ-b^2*x^2*ŷ^2;
f15 = (a*b̂*x̂*y-â*b*x*ŷ)*(a*b*x̂-a*b*ŷ+a*b̂*y-a*x̂*y-â*b*x+â*x*y+b*x*ŷ-b̂*x*y);
eqs = [f1*p[i]^3 * p̂[i]^3 + f2*p[i]^3 * p̂[i]^2 + f3*p[i]^2 * p̂[i]^3 + f4 * p[i]^3 * p̂[i] +  f5*p[i] * p̂[i]^3 + f6 * p[i]^3 + f7 * p̂[i]^3 + f8 * p[i]^2 * p̂[i]^2 + f9 * p[i]^2 * p̂[i] + f10* p[i] * p̂[i]^2 + f11 * p[i]^2 +  f12 * p̂[i]^2 + f13 * p[i] * p̂[i] + f14 * p[i] + f15 *p̂[i] for i=1:8]
F = System(eqs,variables=[x,x̂,ŷ,y,a,â,b,b̂], parameters=vcat(p,p̂))

grading = scaling_symmetries(F)

xp0 = find_start_pair(F)

F = run_monodromy(F, xp0)
symmetries = compute_symmetries!(F, degree=2, param_dep=false)


# Another Formulation
using DecomposingPolynomialSystems, HomotopyContinuation
@var x a y b x̄ ā ȳ b̄
@var γ[1:8] γ̄[1:8]
@var δ[1:8] δ̄[1:8]

D1 = [(ā * x - δ̄[i] * x) * γ[i] + (a * x̄ - δ[i] * x̄) *  γ̄[i] +
      (ā - x̄) * δ[i] + (a - x) * δ̄[i] - δ[i] * δ̄[i] for i in 1:8]
D2 = [(b̄ * y - δ̄[i] * y) * γ[i] + (b * ȳ - δ[i] * ȳ) * γ̄[i] +
      (b̄ - ȳ) * δ[i] + (b - y) * δ̄[i] - δ[i] * δ̄[i] for i in 1:8]
D3 = [γ[i] * γ̄[i] + γ[i] + γ̄[i] for i in 1:8]
eqs = [D1; D2; D3]

F = System(eqs; variables=vcat([x, a, y, b, x̄, ā, ȳ, b̄], γ, γ̄), parameters=vcat(δ, δ̄))
grading = scaling_symmetries(F)

vars = vcat(variables(F), parameters(F))
degree = 2
MDs = multidegrees_up_to_total_degree(length(vars), degree)
classes = partition_multidegrees(MDs, grading)
max(length.(collect(values(classes)))...)

xp0 = HomotopyContinuation.find_start_pair(F)
F = run_monodromy(F, xp0, max_loops_no_progress=1)
F.symmetry_permutations
symmetries = compute_symmetries_graded!(F, grading, MDs, classes)