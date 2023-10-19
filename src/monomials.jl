export Monomial,
    MonomialVector,
    monomials,
    to_expressions,
    to_classes,
    only_param_dep

Grading = Vector{Tuple{Int, Matrix{Int}}}

# TODO: remove?
# TODO: extend Number or nothing at all?
struct Monomial{T<:Integer} <: Number
    md::Vector{T}
    vars::Vector{Variable}

    function Monomial{T}(md, vars) where {T<:Integer}
        # TODO: check if lengths are equal, if mds contains numbers >= 0
        return new(md, vars)
    end
end

mutable struct MonomialVector{T<:Integer} # <: AbstractVector{Monomial}
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

# TODO
function Base.show(io::IO, mons::MonomialVector)
    print(io, "MonomialVector of length $(length(mons.mds))")
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

function monomials(
    vars::Vector{Variable},
    d::Integer;
    homogeneous::Bool=false
)
    homogeneous && return MonomialVector(multidegrees_homogeneous(length(vars), d), vars)
    return MonomialVector(multidegrees_affine(length(vars), d), vars)
end

function to_expressions(mons::MonomialVector)
    nonzero_ids = [findall(!iszero, md) for md in mons.mds]
    return [prod(mons.vars[ids].^md[ids]) for (md, ids) in zip(mds, nonzero_ids)]
end

function Base.gcd(mons::MonomialVector)
    return Monomial(vec(minimum(hcat(mds...); dims=2)), mons.vars)
end

function HC.degree(md::Vector{<:Integer}, grading::Grading)
    return vcat([modV(Uᵢ*md, sᵢ) for (sᵢ, Uᵢ) in grading]...)
end

HC.degree(mon::Monomial, grading::Grading) = degree(mon.md, grading)

# TODO: Can we save multidegrees immediately?
function to_classes(mons::MonomialVector, grading::Grading)
    classes = Dict{Vector{Int}, Vector{Int}}()
    for (i, md) in enumerate(mons.mds)
        deg = HC.degree(md, grading)
        if isnothing(get(classes, deg, nothing)) # the key doesn't exist
            classes[deg] = [i]
        else
            push!(classes[deg], i)
        end
    end
    return classes
end

only_param_dep(md::Vector{<:Integer}, n_unknowns::Integer) = iszero(md[1:n_unknowns])
only_param_dep(mon::Monomial, n_unknowns::Integer) = only_param_dep(mon.md, n_unknowns)
only_param_dep(
    mds::Vector{Vector{<:Integer}},
    n_unknowns::Integer
) = all([only_param_dep(md, n_unknowns) for md in mds])
only_param_dep(mons::MonomialVector, n_unknowns::Integer) = only_param_dep(mons.mds, n_unknowns)
