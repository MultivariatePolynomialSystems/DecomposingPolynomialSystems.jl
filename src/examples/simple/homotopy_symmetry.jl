using DecomposingPolynomialSystems, HomotopyContinuation

@var x a b t
eqs = [(1-t)*(x^2+a) + t*(x^2+a*x+b)]
F = System(eqs; variables=[x], parameters=[a,b,t])

F = run_monodromy(F)
sample_system!(F, 20)

symmetries = compute_symmetries!(F, degree=2, param_dep=true)
