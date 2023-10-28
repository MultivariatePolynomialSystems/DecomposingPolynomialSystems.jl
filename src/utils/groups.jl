export group_structure, to_group, to_permutations

using GAP: GapObj, julia_to_gap, gap_to_julia, Globals
Gl = Globals
to_gap = julia_to_gap
to_julia = gap_to_julia

group_structure(G::GapObj) = String(Gl.StructureDescription(G))

function group_structure(perms::Vector{Vector{Int}})
    return to_julia(Gl.StructureDescription(to_group(perms)))
end

function to_group(perms::Vector{Vector{Int}})
    if length(perms) == 0
        return GAP.evalstr( "Group(())" )  # TODO: raise error?
    end
    Sym = Gl.SymmetricGroup(length(perms[1]))
    gap_gens = [Gl.PermList(to_gap(perm)) for perm in perms]
    Gal = to_gap(gap_gens)
    return Gl.Subgroup(Sym,Gal)
end

function perm_to_list(perm::GapObj, n::Int)
    return [Int(Gl.OnPoints(i, perm)) for i in 1:n]
end

function to_permutations(G::GapObj)
    elems_gap = Gl.Elements(G)
    n = Gl.LargestMovedPoint(Gl.Parent(G))
    return [perm_to_list(elem, n) for elem in elems_gap]
end

centralizer(G::GapObj) = Gl.Centralizer(Gl.Parent(G), G)

function centralizer(perms::Vector{Vector{Int}})
    if length(perms) == 0
        return GAP.evalstr( "Group(())" )
    end
    Sym = Gl.SymmetricGroup(length(perms[1]))
    cents = [Gl.Centralizer(Sym, Gl.PermList(to_gap(perms[i]))) for i in eachindex(perms)]
    return Gl.Intersection(cents...)
end

function block_partition(G::GapObj)
    n = Gl.LargestMovedPoint(G)
    return Vector{Vector{Int}}(to_julia(Gl.Blocks(G, to_gap(Vector(1:n)))))
end

function all_block_partitions(G::GapObj)
    all_blocks = to_julia(Vector{Vector{Int}}, Gl.AllBlocks(G))
    n = Gl.LargestMovedPoint(G)
    block_partitions = Vector{Vector{Vector{Int}}}([])
    for block in all_blocks
        block_partition = Gl.Blocks(G, to_gap(Vector{Int}(1:n)), to_gap(block))
        push!(block_partitions, to_julia(Vector{Vector{Int}}, block_partition))
    end
    return block_partitions
end

function action_on_blocks(G::GapObj, block_partition::Vector{Vector{Int}})
    blocks_gap = to_gap([to_gap(block) for block in block_partition])
    return Gl.Action(G, blocks_gap, Gl.OnSets)
end

function action_on_block(
    G::GapObj,
    block_partition::Vector{Vector{Int}},
    block_id::Int
)
    n_blocks = length(block_partition)
    B = block_partition[block_id]
    block_size = length(B)
    elems = Set(vcat(block_partition...))
    remaining = collect(setdiff(elems, Set(B)))
    new_blocks = to_gap([to_gap(sort(vcat([B[j]], remaining))) for j in 1:block_size])
    S = Gl.Stabilizer(G, to_gap(B), Gl.OnSets)
    return Gl.Action(S, new_blocks, Gl.OnSets)
end

function action_on_given_blocks(
    G::GapObj,
    block_partition::Vector{Vector{Int}},
    block_ids::Vector{Int}
)
    Bs = sort(vcat(block_partition[block_ids]...))
    set_size = length(Bs)
    elems = Set(vcat(block_partition...))
    remaining = collect(setdiff(elems, Set(Bs)))
    new_blocks = to_gap([to_gap(sort(vcat([Bs[j]], remaining))) for j in 1:set_size])
    S = Gl.Stabilizer(G, to_gap(Bs), Gl.OnSets)
    return Gl.Action(S, new_blocks, Gl.OnSets)
end

# Intersection of all block stabilizers
function kernel_of_action_on_blocks(G::GapObj, block_partition::Vector{Vector{Int}})
    B1 = block_partition[1]
    K = Gl.Stabilizer(G, to_gap(B1), Gl.OnSets)
    n_blocks = length(block_partition)
    for i = 2:n_blocks
        Bi = block_partition[i]
        Si = Gl.Stabilizer(G, to_gap(Bi), Gl.OnSets)
        K = Gl.Intersection(K, Si)
        println("i = ", i)
    end
    return K
end
