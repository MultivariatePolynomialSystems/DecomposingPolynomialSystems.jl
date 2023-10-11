using DecomposingPolynomialSystems

@var x[1:2] p[1:2]
F = System([x[1]^2 - x[2]^2 - p[1], 2*x[1]*x[2] - p[2]]; variables=x, parameters=p)

F = run_monodromy(F)

deck = symmetries_fixing_parameters!(F; degree_bound=1, param_dep=false)
deck[2]