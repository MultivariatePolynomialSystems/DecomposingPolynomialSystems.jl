export Automorphism,
    AutomorphismGroup

MiExpression = Union{Missing, Expression}

struct Automorphism
    exprs::Vector{MiExpression}
    unknowns::Vector{Variable}
    parameters::Vector{Variable}

    function Automorphism(exprs, unknowns, parameters)
        # TODO: verify args
        return new(exprs, unknowns, parameters)
    end
end

Base.getindex(dt::Automorphism, inds...) = getindex(dt.exprs, inds...)

function Base.show(io::IO, dt::Automorphism)
    println(
        io,
        "Automorphism: acts on $(phrase(length(dt.unknowns), "unknown")),",
        " fixes $(phrase(length(dt.parameters), "parameter"))",
    )
    println(io, " action:")
    for i in 1:length(dt.exprs)
        print(io, "  ", dt.unknowns[i], " ↦ ", dt.exprs[i])
        i < length(dt.exprs) && print(io, "\n")
    end
end

"""
    AutomorphismGroup

A `AutomorphismGroup` is the result of automorphisms computation.
"""
struct AutomorphismGroup
    maps::Vector{Automorphism}
    group::GapObj
    F::SampledParametricSystem
end

function AutomorphismGroup(F::SampledParametricSystem)
    symmetries = _init_automorphisms(length(aut_permutations(F)), unknowns(F))
    return AutomorphismGroup(symmetries, F)
end

function AutomorphismGroup(
    symmetries::Vector{Vector{MiExpression}},
    F::SampledParametricSystem
)
    action = [Automorphism(symmetry, unknowns(F), parameters(F)) for symmetry in symmetries]
    return AutomorphismGroup(action, to_group(aut_permutations(F)), F)
end

function Base.show(io::IO, deck::AutomorphismGroup)
    println(io, "AutomorphismGroup of order $(length(deck.maps))")
    println(io, " structure: ", order(deck.group) == 1 ? "trivial" : group_structure(deck.group))
    print(io, " action:")
    for i in eachindex(deck.maps)
        println(io, "\n  ", to_ordinal(i), " map:")
        for (j, var) in enumerate(unknowns(deck.F))  # action on parameters is trivial, don't show it
            print(io, "   ", var, " ↦ ", deck.maps[i][j])
            j < length(unknowns(deck.F)) && print(io, "\n")
        end
    end
end

Base.getindex(deck::AutomorphismGroup, inds...) = getindex(deck.maps, inds...)
