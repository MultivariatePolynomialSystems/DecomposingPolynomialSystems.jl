export Monomial,
    MonomialVector,
    to_expressions,
    to_classes,
    only_param_dep,
    n_only_param_dep

# TODO: remove?
# TODO: extend Number or nothing at all?
struct Monomial{T<:Integer} #<: Number
    md::Vector{T}
    vars::Vector{Variable}

    function Monomial{T}(md, vars) where {T<:Integer}
        # TODO: check if lengths are equal, if mds contains numbers >= 0
        return new(md, vars)
    end
end

Monomial(md::Vector{T}, vars::Vector{Variable}) where {T<:Integer} = Monomial{T}(md, vars)
function Monomial{T}(var::Variable, vars::Vector{Variable}) where {T<:Integer}
    md = zeros(T, length(vars))
    md[findfirst(x->x==var, vars)] = 1
    return Monomial{T}(md, vars)
end

function Monomial{T}(
    mon::Expression,
    vars::Vector{Variable}
) where {T<:Integer}
    es, cs = exponents_coefficients(mon, vars)
    if length(cs) > 1 || !isone(cs[1])
        throw(ArgumentError("Input expression is not a monomial"))
    else
        return Monomial{T}(es[:,1], vars)
    end
end

Base.isone(mon::Monomial) = iszero(mon.md)
Base.convert(::Type{Expression}, mon::Monomial) = prod(mon.vars.^mon.md)
Base.:(==)(m₁::Monomial, m₂::Monomial) = Expression(m₁) == Expression(m₂)
Base.show(io::IO, mon::Monomial) = show(io, Expression(mon))

mutable struct MonomialVector{T<:Integer} # <: AbstractVector{Expression}
    mds::Vector{Vector{T}}
    vars::Vector{Variable}

    function MonomialVector{T}(mds, vars) where {T<:Integer}
        # TODO: check if lengths are equal, if mds contains numbers >= 0
        return new(mds, vars)
    end
end

MonomialVector(
    mds::Vector{Vector{T}},
    vars::Vector{Variable}
) where {T<:Integer} = MonomialVector{T}(mds, vars)

Base.length(mons::MonomialVector) = length(mons.mds)
Base.getindex(mons::MonomialVector, i::Integer) = Monomial(mons.mds[i], mons.vars)
Base.getindex(
    mons::MonomialVector,
    inds...
) = MonomialVector(getindex(mons.mds, inds...), mons.vars)
# Base.findfirst(f::Function, mons::MonomialVector) = findfirst(f, mons.mds)
Base.findfirst(m::Monomial, mons::MonomialVector) = findfirst(x->x==m.md, mons.mds)
function Base.push!(mons::MonomialVector{T}, md::Vector{T}) where {T<:Integer}
    push!(mons.mds, md)
end

function Base.vcat(monVs::MonomialVector...)
    # TODO: check if mons.vars are all equal
    MonomialVector(vcat([mons.mds for mons in monVs]...), monVs[1].vars)
end

# TODO
function Base.show(io::IO, mons::MonomialVector)
    # println(io, "$(length(mons.mds))-element $(typeof(mons))")
    print(io, "[", join(to_expressions(mons), ", "), "]")
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
    vars::Vector{Variable};
    degree::Integer,
    homogeneous::Bool=false
) where {T<:Integer}

    degree = convert(T, degree)
    if homogeneous
        return MonomialVector(multidegrees_homogeneous(length(vars), degree), vars)
    end
    return MonomialVector(multidegrees_affine(length(vars), degree), vars)
end

function to_expressions(mons::MonomialVector)
    nonzero_ids = [findall(!iszero, md) for md in mons.mds]
    return [prod(mons.vars[ids].^md[ids]) for (md, ids) in zip(mons.mds, nonzero_ids)]
end

function Base.gcd(mons::MonomialVector)
    return Monomial(vec(minimum(hcat(mons.mds...); dims=2)), mons.vars)
end

only_param_dep(md::Vector{<:Integer}, unknowns_ids::Vector{Int}) = iszero(md[unknowns_ids])
only_param_dep(mon::Monomial, unknowns_ids::Vector{Int}) = only_param_dep(mon.md, unknowns_ids)
only_param_dep(
    mds::Vector{Vector{T}},
    unknowns_ids::Vector{Int}
) where {T<:Integer} = all([only_param_dep(md, unknowns_ids) for md in mds])
only_param_dep(
    mons::MonomialVector,
    unknowns_ids::Vector{Int}
) = only_param_dep(mons.mds, unknowns_ids)

n_only_param_dep(
    mons::MonomialVector,
    unknown_ids::Vector{Int}
) = sum([only_param_dep(md, unknown_ids) for md in mons.mds])