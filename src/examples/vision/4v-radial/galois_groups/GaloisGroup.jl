include("../../../src/utils/utils.jl")
include("permutations.jl")

ps = [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13]

G = permutations_to_group(ps)
block_partitions = all_block_partitions(G)

# Galois group of calibrated formulation (3584 sols)
G1 = G
B1 = block_partition(G1)
# GAP.Globals.Order(G1)
# oG1 = factorial(big(28))*2^27*(factorial(big(4)))^28*2^56*big(2)^364

G2 = action_on_blocks(G1, B1)
B2 = block_partition(G2)

G3 = action_on_blocks(G2, B2)
B3 = block_partition(G3)

G4 = action_on_blocks(G3, B3)
B4 = block_partition(G4)

# Galois group of metric upgrade (224 sols)
G5 = action_on_blocks(G4, B4)
B5 = block_partition(G5)
GAP.Globals.Order(G5)
# oG5 = factorial(big(28))*2^27*(factorial(big(4)))^28*2^56
#
# AonB_G5 = action_on_block(G5, B5, 1)
# GAP.Globals.Order(AonB_G5)
#
# A224 = GAP.Globals.AlternatingGroup(224)
# GAP.Globals.IsSubgroup(A224, G5)

A5 = GAP.Globals.ActionHomomorphism(G5, GAP.julia_to_gap([GAP.julia_to_gap(B5[i]) for i=1:length(B5)]), GAP.Globals.OnSets)
K = GAP.Globals.Kernel(A5)
# K1 of order 4!^28 --> S4 acting on 28 pairs of blocks, on 2 blocks from the same pair acts in the same way
# K2 of order 2^56 --> transposition in every block or
# K = K1 x K2
Kgen = GAP.Globals.GeneratorsOfGroup(K)

function find_block(block_partition::Vector{Vector{Int}}, elem_id::Int)
    for i=1:length(block_partition)
        if elem_id in block_partition[i]
            return i
        end
    end
    return 0
end

Kgen[1]
Kgen[2]
Kgen[3]
Kgen[4]
Kgen[5]
Kgen[6]
Kgen[7]
Kgen[8]
Kgen[9]
Kgen[10]
Kgen[11]
Kgen[12]
Kgen[13]
Kgen[14]
Kgen[15]
Kgen[16]
Kgen[17]
Kgen[18]
Kgen[19]
Kgen[20]



find_block(B5, 108)
find_block(B5, 160)
GAP.julia_to_gap(B5[29])
GAP.julia_to_gap(B5[56])
# S4 seems to act on the following pairs of blocks
transpositions_and_cycle3_block_pairs = [[42, 54], [33, 53], [51, 56], []]
cycle3_block_pairs = [[20, 50], [22, 33], [29, 56], []]
# involutions inside blocks: (108,157)(156,178), (108,178)(156,157)


# Kg = GAP.Globals.MinimalGeneratingSet(K)

find_block(B5, 136)
find_block(B5, 202)

[[140, 193], [158, 215]] # B54 & B42
[[156, 178], [164, 182]] # B55 & B18
[[120, 142], [123, 176]] # B12 & B41
[[124, 187], [145, 217]]
[[147, 175], [109, 137]]
[[186, 206], [99, 104]]
[[159, 222], [170, 218]]
[[138, 207], [171, 190]]
[[136, 151], [202, 208]]









# Galois group of uncalibrated formulation (56 sols)
G6 = action_on_blocks(G5, B5)
B6 = block_partition(G6)
GAP.Globals.Order(G6)
oG6 = factorial(big(28))*2^27

A56 = GAP.Globals.AlternatingGroup(56)
GAP.Globals.IsSubgroup(A56, G6)

# Galois group of tensor formulation (28 sols)
G7 = action_on_blocks(G6, B6)
B7 = block_partition(G7)
GAP.Globals.Order(G7)
oG7 = factorial(big(28))
