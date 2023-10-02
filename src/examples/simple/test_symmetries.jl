using DecomposingPolynomialSystems, HomotopyContinuation

@var x[1:2] p[1:2]
F = System([x[1]^2 - x[2]^2 - p[1], 2*x[1]*x[2] - p[2]]; variables=x, parameters=p)

F = run_monodromy(F)
DeckTransformationGroup(F)

symmetries_fixing_parameters_dense!(F; degree_bound=1, param_dep=false)

scalings = scaling_symmetries(F, x)
scalings.grading[1]
scalings.grading[2]

mds = multidegrees_up_to_total_degree(2, 1)
classes = partition_multidegrees(mds, scalings.grading)

F.symmetry_permutations
permutations_to_group(F.symmetry_permutations)

d = Dict(zip([1,2,3], [4,5,6]))
for (a, d[a]) in d
    print(a)
end

