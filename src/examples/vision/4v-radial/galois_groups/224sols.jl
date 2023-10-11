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

# Galois group of metric upgrade (224 sols)
G5 = action_on_blocks(G4, B4)
B5 = block_partition(G5)

# Action on 56 blocks of size 4
A = GAP.Globals.ActionHomomorphism(G5, GAP.julia_to_gap([GAP.julia_to_gap(B5[i]) for i=1:length(B5)]), GAP.Globals.OnSets)
K = GAP.Globals.Kernel(A)

A = action_on_given_blocks(K, B5, [1,2])
GAP.Globals.StructureDescription(A)

block_pairs = [[1,39], [2,44], [3,43], [4,49], [5,34],
               [6, 21], [7,31], [8,45], [9,47], [10,27], [11,35], [12,41],
               [13,52], [14,40], [15,37], [16,26], [17,30],
               [18,55], [19,36], [20,50], [22,33], [23,48], [24,51],
               [25,38], [28,46], [29,56], [32,53], [42,54]]
for pair in block_pairs
    A = action_on_given_blocks(K, B5, pair)
    println(GAP.Globals.StructureDescription(A))
end

A = action_on_given_blocks(K, B5, block_pairs[1])
GAP.Globals.StructureDescription(A)
GAP.Globals.Order(A) # 96 = 4 * 4!
GAP.Globals.Order(K)
big(96)^28

# Proving that A is isomorphic to (S2 x S2) \rtimes S4
Ns = [GAP.Globals.SymmetricGroup(4),
      GAP.Globals.SymmetricGroup(4),
      GAP.Globals.CyclicGroup(4),
      GAP.Globals.DirectProduct(GAP.Globals.CyclicGroup(2), GAP.Globals.CyclicGroup(2))]
Hs = [GAP.Globals.DirectProduct(GAP.Globals.SymmetricGroup(2), GAP.Globals.SymmetricGroup(2)),
      GAP.Globals.CyclicGroup(4),
      GAP.Globals.SymmetricGroup(4),
      GAP.Globals.SymmetricGroup(4)]

i = 4
N, H = Ns[i], Hs[i]
AutN = GAP.Globals.AutomorphismGroup(N)
homs = GAP.Globals.AllHomomorphisms(H, AutN)
for j=1:GAP.Globals.Length(homs)
    NsemiH = GAP.Globals.SemidirectProduct(H, homs[j], N)
    println("Are isomorphic: ", GAP.Globals.IsomorphismGroups(NsemiH, A))
end


# Hypothesis: K = K1 x ... x K28, where Kj is the stabilizer of the j-th pair of blocks and all the remaining elements

# Action on the set consisting of 224-8+1=217 elements: one pair of blocks and the remaining elements
for j=1:28
    S1 = sort(vcat(B5[block_pairs[j]]...))

    r = [n for n=1:28 if n != j]
    L = vcat(B5[vcat(block_pairs[r]...)]...)
    S2 = [[n] for n in L]
    S = vcat([S1], S2)

    A = GAP.Globals.ActionHomomorphism(K, GAP.julia_to_gap([GAP.julia_to_gap(S[i]) for i=1:length(S)]), GAP.Globals.OnSets)

    # Elements in K which stabilize the first pair of blocks and every of the remaining elements
    Kj = GAP.Globals.Kernel(A)
    println("Order = ", GAP.Globals.Order(Kj))
    println(GAP.Globals.SmallGeneratingSet(Kj))
end

# Let H,K \subset G. Then
# G \cong H x K if and only if
# 1) G = HK
# 2) H \cap K = {id}
# 3) H,K < G (normal subgroups), or hk = kh for all h \in H, k \in K

# 1) Obviously, K1*...*K28 < K. Also, K < K1*...*K28 since |K1*...*K28| = |K|
# 2) Ki \cap Kj = {id} for i \neq j, since every permutation from the intersection
# fixes the elements S1 = [224]\{Bi1,Bi2} (i.e. it lies in Ki) and the elements
# S2 = [224]\{Bj1,Bj2} (i.e. it lies in Kj). Since i \neq j, then S1 \cup S2 = [224],
# i.e. p is the identity permutation.
# 3) To show that Kj is a normal subgroup of K we need to prove that
# pKj = Kjp for every p \in K. If p \in Kj, it is trivial. If p \in Ki, i \neq j,
# then take any a \in Kj. If p = p1p2...pn and a = a1a2...ak are representations
# of p and a as products of disjoint cycles, then obiously cycles
# p1,...,pn,a1,...,ak are disjoint as well, since they permute elements in
# disjoint i-th and j-th pairs of blocks.
