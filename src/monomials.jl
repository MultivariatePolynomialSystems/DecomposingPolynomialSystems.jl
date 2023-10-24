export Monomial,
    MonomialVector,
    monomials,
    to_expressions,
    to_classes,
    only_param_dep

# TODO: remove?
# TODO: extend Number or nothing at all?
struct Monomial{T<:Integer} #<: Number
    vars::Vector{Variable}
    md::Vector{T}

    function Monomial{T}(vars, md) where {T<:Integer}
        # TODO: check if lengths are equal, if mds contains numbers >= 0
        return new(vars, md)
    end
end

Monomial(vars::Vector{Variable}, md::Vector{T}) where {T<:Integer} = Monomial{T}(vars, md)
Base.isone(mon::Monomial) = iszero(mon.md)
Base.convert(::Type{Expression}, mon::Monomial) = prod(mon.vars.^mon.md)
Base.show(io::IO, mon::Monomial) = show(io, Expression(mon))

mutable struct MonomialVector{T<:Integer} # <: AbstractVector{Expression}
    vars::Vector{Variable}
    mds::Vector{Vector{T}}

    function MonomialVector{T}(vars, mds) where {T<:Integer}
        # TODO: check if lengths are equal, if mds contains numbers >= 0
        return new(vars, mds)
    end
end

MonomialVector(
    vars::Vector{Variable},
    mds::Vector{Vector{T}}
) where {T<:Integer} = MonomialVector{T}(vars, mds)

Base.length(mons::MonomialVector) = length(mons.mds)
Base.getindex(mons::MonomialVector, i::Integer) = Monomial(mons.mds[i], mons.vars)
Base.getindex(
    mons::MonomialVector,
    inds...
) = MonomialVector(mons.vars, getindex(mons.mds, inds...))

function Base.vcat(monVs::MonomialVector...)
    # TODO: check if mons.vars are all equal
    MonomialVector(monVs[1].vars, vcat([mons.mds for mons in monVs]...))
end

# TODO
function Base.show(io::IO, mons::MonomialVector)
    print(io, "$(typeof(mons)) of length $(length(mons.mds))")
end

function multidegrees_homogeneous(md::Vector{T}, n::Integer, d::T) where {T<:Integer}
    n == 1 && return [vcat(md, d)]
    return vcat([multidegrees_homogeneous(vcat(md, convert(T, i)), n-1, d-i) for i::T in 0:d]...)
end

function multidegrees_homogeneous(n::Integer, d::T) where {T<:Integer}
    return multidegrees_homogeneous(Vector{T}([]), n, d)
end

function multidegrees_affine(n::Integer, d::T) where {T<:Integer}
    return vcat([multidegrees_homogeneous(n, i) for i::T in 0:d]...)
end

function MonomialVector{T}(
    vars::Vector{Variable},
    d::Integer;
    homogeneous::Bool=false
) where {T<:Integer}

    d = convert(T, d)
    if homogeneous
        return MonomialVector(vars, multidegrees_homogeneous(length(vars), d))
    end
    return MonomialVector(vars, multidegrees_affine(length(vars), d))
end

function to_expressions(mons::MonomialVector)
    nonzero_ids = [findall(!iszero, md) for md in mons.mds]
    return [prod(mons.vars[ids].^md[ids]) for (md, ids) in zip(mons.mds, nonzero_ids)]
end

function Base.gcd(mons::MonomialVector)
    return Monomial(vec(minimum(hcat(mons.mds...); dims=2)), mons.vars)
end

# Methods below suppose that vars = [unknowns, parameters]
only_param_dep(md::Vector{<:Integer}, n_unknowns::Integer) = iszero(md[1:n_unknowns])
only_param_dep(mon::Monomial, n_unknowns::Integer) = only_param_dep(mon.md, n_unknowns)
only_param_dep(
    mds::Vector{Vector{T}},
    n_unknowns::Integer
) where {T<:Integer} = all([only_param_dep(md, n_unknowns) for md in mds])
only_param_dep(mons::MonomialVector, n_unknowns::Integer) = only_param_dep(mons.mds, n_unknowns)
