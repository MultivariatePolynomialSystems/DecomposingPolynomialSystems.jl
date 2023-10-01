using DecomposingPolynomialSystems, HomotopyContinuation

@var x[1:2] p[1:2]
F = System([x[1]^2 - x[2]^2 - p[1], 2*x[1]*x[2] - p[2]]; variables=x, parameters=p)

F = run_monodromy(F)

deck = symmetries_fixing_parameters!(F, degree_bound=1, param_dep=false)
deck.action

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

using DecomposingPolynomialSystems, GAP

perms = [[2, 1, 3, 4]]
g = permutations_to_group(perms)
group_to_permutations(g)
GAP.gap_to_julia(group_structure(g))

GAP.Globals.Parent(g)
elems = GAP.Globals.Elements(g)

i = GAP.Globals.OnPoints(3, elems[2])
GAP.Globals.ListPerm(elems[2])

GAP.Globals.ListPerm(elems[2])
group_to_permutations(g)

GAP.Globals.PermList(GAP.julia_to_gap([1,2,4,3]))