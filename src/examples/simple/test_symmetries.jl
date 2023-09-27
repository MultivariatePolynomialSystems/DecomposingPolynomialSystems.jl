using DecomposingPolynomialSystems, HomotopyContinuation

@var x[1:2] p[1:2]
F = System([x[1]^2 - x[2]^2 - p[1], 2*x[1]*x[2] - p[2]]; variables=x, parameters=p)

deck = symmetries_fixing_parameters(F, degree_bound=1, param_dep=false)
deck.action

scalings = scaling_symmetries(F, x)
scalings.grading[1]
scalings.grading[2]

mds = multidegrees_up_to_total_degree(2, 1)
classes = partition_multidegrees(mds, scalings.grading)