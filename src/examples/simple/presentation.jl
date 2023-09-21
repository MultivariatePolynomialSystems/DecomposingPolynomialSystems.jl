using DecomposingPolynomialSystems

@var x a b
F = System([a*x^4 + b*x^3 + x^2 + b*x + a]; parameters = [a, b])

# Monodromy group
F = run_monodromy(F)
G = permutations_to_group(F.monodromy_permutations)
using GAP
GAP.Globals.StructureDescription(G)
GAP.Globals.Order(G)

# Symmetry group
A = centralizer(F.monodromy_permutations)
GAP.Globals.StructureDescription(A)

# Symmetries
symmetries = compute_symmetries!(F, degree=1)