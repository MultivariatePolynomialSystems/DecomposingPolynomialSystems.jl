include("../../../../src/utils/utils.jl")
include("permutations.jl")

ps = [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13]

G = permutations_to_group(ps)

centralizer(G)

block_partitions = all_block_partitions(G)

B = block_partitions[6]

# Action on 224 blocks of size 16
A = GAP.Globals.ActionHomomorphism(G, GAP.julia_to_gap([GAP.julia_to_gap(B[i]) for i=1:length(B)]), GAP.Globals.OnSets)
K = GAP.Globals.Kernel(A)


A = action_on_given_blocks(K, B, [1,161,190,207])
GAP.Globals.StructureDescription(A)

# 364 = 2*2*91
GAP.Globals.StructureDescription(K)
gens = GAP.Globals.GeneratorsOfGroup(K)

function find_block(block_partition::Vector{Vector{Int}}, elem_id::Int)
    for i=1:length(block_partition)
        if elem_id in block_partition[i]
            return i
        end
    end
    return 0
end

MPs = GAP.gap_to_julia(GAP.Globals.MovedPoints(gens[10]))
Set([find_block(B, j) for j in MPs])
