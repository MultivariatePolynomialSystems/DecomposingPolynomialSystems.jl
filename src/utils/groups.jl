export group_structure, permutations_to_group

function filter_permutations(perms::Matrix{Int})::Vector{Vector{Int}}
    nsols = length(perms[:,1])
    return filter!(x->!(0 in x) & (length(unique(x)) == nsols), [perms[:,i] for i in 1:size(perms, 2)])
end

function group_structure(G::GapObj)
    return GAP.Globals.StructureDescription(G)
end

function permutations_to_group(perms::Vector{Vector{Int}})::GapObj
    if length(perms) == 0
        return GAP.evalstr( "Group(())" )
    end
    Sym = GAP.Globals.SymmetricGroup(length(perms[1]))
    gap_gens=[]
    for i in eachindex(perms)
        push!(gap_gens, GAP.Globals.PermList(GAP.julia_to_gap(perms[i])));
    end
    Gal = GAP.julia_to_gap(gap_gens)
    return GAP.Globals.Subgroup(Sym,Gal)
end

function block_partition(G::GapObj)
    n = GAP.Globals.LargestMovedPoint(G) # TODO
    return Vector{Vector{Int64}}(GAP.gap_to_julia(GAP.Globals.Blocks(G, GAP.julia_to_gap(Vector{Int64}(1:n)))))
end

function all_block_partitions(G::GapObj)::Vector{Vector{Vector{Int}}}
    all_blocks = GAP.gap_to_julia(Vector{Vector{Int64}}, GAP.Globals.AllBlocks(G))
    n = GAP.Globals.LargestMovedPoint(G) # TODO
    block_partitions = []
    for block in all_blocks
        block_partition = GAP.Globals.Blocks(G, GAP.julia_to_gap(Vector{Int64}(1:n)), GAP.julia_to_gap(block))
        push!(block_partitions, GAP.gap_to_julia(Vector{Vector{Int64}}, block_partition))
    end
    return Vector{Vector{Vector{Int64}}}(block_partitions)
end

function centralizer(G::GapObj)::GapObj
    return GAP.Globals.Centralizer(GAP.Globals.Parent(G), G)
end

function centralizer(perms::Vector{Vector{Int}})::GapObj
    if length(perms) == 0
        return GAP.evalstr( "Group(())" )
    end
    Sym = GAP.Globals.SymmetricGroup(length(perms[1]))
    cents=[]
    for i in eachindex(perms)
        push!(cents, GAP.Globals.Centralizer(Sym, GAP.Globals.PermList(GAP.julia_to_gap(perms[i]))));
    end
    return GAP.Globals.Intersection(cents...)
end

function group_to_permutations(G::GapObj)::Vector{Vector{Int}}
    elems_gap = GAP.Globals.Elements(G)
    elems_julia = []
    for i = 1:GAP.Globals.Length(elems_gap)
        push!(elems_julia, GAP.gap_to_julia(Vector{Int64}, GAP.Globals.ListPerm(elems_gap[i])))
    end
    return elems_julia
end

function action_on_blocks(G::GapObj, block_partition::Vector{Vector{Int}})::GapObj
    blocks_gap = GAP.julia_to_gap([GAP.julia_to_gap(block) for block in block_partition])
    return GAP.Globals.Action(G, blocks_gap, GAP.Globals.OnSets)
end

function action_on_block(G::GapObj, block_partition::Vector{Vector{Int}}, block_id::Int)::GapObj
    n_blocks = length(block_partition)
    B = block_partition[block_id]
    block_size = length(B)
    elems = Set(vcat(block_partition...))
    remaining = collect(setdiff(elems, Set(B)))
    new_blocks = GAP.julia_to_gap([GAP.julia_to_gap(sort(vcat([B[j]], remaining))) for j in 1:block_size])
    S = GAP.Globals.Stabilizer(G, GAP.julia_to_gap(B), GAP.Globals.OnSets)
    return GAP.Globals.Action(S, new_blocks, GAP.Globals.OnSets)
end

function action_on_given_blocks(G::GapObj, block_partition::Vector{Vector{Int}}, block_ids::Vector{Int})::GapObj
    Bs = sort(vcat(block_partition[block_ids]...))
    set_size = length(Bs)
    elems = Set(vcat(block_partition...))
    remaining = collect(setdiff(elems, Set(Bs)))
    new_blocks = GAP.julia_to_gap([GAP.julia_to_gap(sort(vcat([Bs[j]], remaining))) for j in 1:set_size])
    S = GAP.Globals.Stabilizer(G, GAP.julia_to_gap(Bs), GAP.Globals.OnSets)
    return GAP.Globals.Action(S, new_blocks, GAP.Globals.OnSets)
end

# Intersection of all block stabilizers
function kernel_of_action_on_blocks(G::GapObj, block_partition::Vector{Vector{Int}})::GapObj
    B1 = block_partition[1]
    K = GAP.Globals.Stabilizer(G, GAP.julia_to_gap(B1), GAP.Globals.OnSets)
    n_blocks = length(block_partition)
    for i = 2:n_blocks
        Bi = block_partition[i]
        Si = GAP.Globals.Stabilizer(G, GAP.julia_to_gap(Bi), GAP.Globals.OnSets)
        K = GAP.Globals.Intersection(K, Si)
        println("i = ", i)
    end
    return K
end
