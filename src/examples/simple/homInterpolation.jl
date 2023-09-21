using DecomposingPolynomialSystems, HomotopyContinuation
@var x a b c
F = System([x^2+b*x+c^2]; variables=[x], parameters=[b,c])

grading = scaling_symmetries(F)
U, s = grading.U, grading.s
U[1], s[1]
grading = Grading([U[1]], [s[1]])

n_vars = 3
degree = 1
MDs = multidegreesUpToTotalDegree(n_vars, degree)
vars = vcat(variables(F), parameters(F))
classes = partitionMultidegrees(MDs, grading, vars)

F = run_monodromy(F)
symmetries = compute_symmetries_graded!(F, grading, classes)



# (x,b,c) --> (-b-x, b, c) --> (b+x,-b,c)

# (x,b,c) --> (-x,-b,c) --> (b+x,-b,c)



using DecomposingPolynomialSystems
@var x a b c
F = System([a*x^2+b*x+c]; variables=[x], parameters=[a,b,c])

compute_symmetries!(F, degree=3) # smth strange is happenning for degree=2





using DecomposingPolynomialSystems, HomotopyContinuation

@var x a b c
F = System([a*x^2+b*x+c]; variables=[x], parameters=[a,b,c])
F = run_monodromy(F)

U, s = snfScalingSymmetries(F.equations)

classes = compute_symmetries_graded!(F, U, s; degree=2)
