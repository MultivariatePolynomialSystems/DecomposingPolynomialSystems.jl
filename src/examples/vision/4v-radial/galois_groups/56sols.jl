include("../../../../src/utils/utils.jl")
include("permutations.jl")

ps = [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13]

# Galois group of calibrated formulation (3584 sols)
G1 = permutations_to_group(ps)
B1 = block_partition(G1)

G2 = action_on_blocks(G1, B1)
B2 = block_partition(G2)

G3 = action_on_blocks(G2, B2)
B3 = block_partition(G3)

G4 = action_on_blocks(G3, B3)
B4 = block_partition(G4)

G5 = action_on_blocks(G4, B4)
B5 = block_partition(G5)

# Galois group of uncalibrated formulation (56 sols)
G6 = action_on_blocks(G5, B5)
B6 = block_partition(G6)

# Action on 28 blocks of size 2
A = GAP.Globals.ActionHomomorphism(G6, GAP.julia_to_gap([GAP.julia_to_gap(B6[i]) for i=1:length(B6)]), GAP.Globals.OnSets)
K = GAP.Globals.Kernel(A)

GAP.Globals.StructureDescription(K)
gens = GAP.Globals.GeneratorsOfGroup(K)

gens[1]

A56 = GAP.Globals.AlternatingGroup(56)
GAP.Globals.IsSubgroup(A56, K)
