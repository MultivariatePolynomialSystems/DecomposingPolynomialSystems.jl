export Monomial,
    MonomialVector,
    monomials

# TODO: remove?
# TODO: extend Number or nothing at all?
struct Monomial{T<:Integer} <: Number
    md::Vector{T}
    vars::Vector{Variable}

    function Monomial{T}(md::Vector{T}, vars::Vector{Variable}) where {T<:Integer}
        # TODO: check if lengths are equal, if mds contains numbers >= 0
        return new(md, vars)
    end
end

struct MonomialVector{T<:Integer} <: AbstractVector{Monomial}
    mds::Vector{Vector{T}}
    vars::Vector{Variable}

    function MonomialVector{T}(mds::Vector{Vector{T}}, vars::Vector{Variable}) where {T<:Integer}
        # TODO: check if lengths are equal, if mds contains numbers >= 0
        return new(mds, vars)
    end
end

Base.length(mons::MonomialVector) = length(mons.mds)
Base.getindex(mons::MonomialVector, i::Integer) = Monomial(mons.mds[i], mons.vars)
Base.getindex(mons::MonomialVector, inds...) = MonomialVector(getindex(mons.mds, inds...), mons.vars)

# TODO
function Base.show(io::IO, mons::MonomialVector)
    println(io, "MonomialVector")
end

function multidegrees_homogeneous(md::Vector{T}, n::Integer, d::T) where {T<:Integer}
    n == 1 && return [vcat(md, d)]
    return vcat([multidegrees_homogeneous(vcat(md, convert(T, i)), n-1, d-i) for i in 0:d]...)
end

function multidegrees_homogeneous(n::Integer, d::T) where {T<:Integer}
    return multidegrees_homogeneous(Vector{T}([]), n, d)
end

function multidegrees_affine(n::Integer, d::T) where {T<:Integer}
    return vcat([multidegrees_homogeneous(n, convert(T, i)) for i in 0:d]...)
end

function monomials(vars::Vector{Variable}, d::Integer; homogeneous::Bool=false)
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

function HomotopyContinuation.degree(md::Multidegree, grading::Grading)
    return vcat([modV(Uᵢ*md, sᵢ) for (sᵢ, Uᵢ) in grading]...)
end

HomotopyContinuation.degree(mon::Monomial, grading::Grading) = degree(mon.md, grading)

# TODO: Can we save multidegrees immediately?
function to_classes(mds::Vector{Multidegree}, grading::Grading)
    classes = Dict{Vector{Int}, Vector{Int}}()
    for (i, md) in enumerate(mds)
        deg = degree_wrt_grading(md, grading)
        if isnothing(get(classes, deg, nothing)) # the key doesn't exist
            classes[deg] = [i]
        else
            push!(classes[deg], i)
        end
    end
    return classes
end

only_param_dep(md::Multidegree, n_unknowns::Integer) = iszero(md[1:n_unknowns])
only_param_dep(mon::Monomial, n_unknowns::Integer) = only_param_dep(mon.md, n_unknowns)
only_param_dep(
    mds::Vector{Multidegree},
    n_unknowns::Integer
) = all([only_param_dep(md, n_unknowns) for md in mds])
only_param_dep(mons::MonomialVector, n_unknowns::Integer) = only_param_dep(mons.mds, n_unknowns)

# supposes each md in mds is a multidegree in both unknowns and parameters
# TODO: The implementation below (with _) is more efficient (approx 2x),
# TODO: since it exploits the sparsity of multidegrees. REMOVE THIS METHOD?
function HomotopyContinuation.evaluate(mons::MonomialVector, samples::VarietySamples)
    solutions = samples.solutions
    parameters = samples.parameters

    n_unknowns, n_sols, n_instances = size(solutions)
    mds = mons.mds
    n_mds = length(mds)

    evaluated_mons = zeros(CC, n_mds, n_sols, n_instances)
    for i in 1:n_instances
        params = view(parameters, :, i)
        params_eval = [prod(params.^md[n_unknowns+1:end]) for md in mds]
        sols = view(solutions, :, :, i)
        for j in 1:n_mds
            evaluated_mons[j, :, i] = vec(prod(sols.^mds[j][1:n_unknowns], dims=1)).*params_eval[j]
        end
    end
    return evaluated_mons
end

# TODO: consider view for slices
function evaluate_monomials_at_samples_(
    mons::MonomialVector,
    samples::VarietySamples;
    sparse::Bool=false
)
    solutions = samples.solutions
    parameters = samples.parameters

    n_unknowns, n_sols, n_instances = size(solutions)
    mds = mons.mds
    n_mds = length(mds)

    nonzero_ids_unknowns = [findall(!iszero, md[1:n_unknowns]) for md in mds]
    nonzero_ids_params = [findall(!iszero, md[n_unknowns+1:end]) for md in mds]

    evaluated_mons = zeros(CC, n_mds, n_sols, n_instances)
    for i in 1:n_instances
        params = view(parameters, :, i)
        sols = view(solutions, :, :, i)
        for (j, md) in enumerate(mds)
            params_part = params[nonzero_ids_params[j]]
            md_params_part = md[n_unknowns+1:end][nonzero_ids_params[j]]
            params_eval = prod(params_part.^md_params_part)
            sols_part = view(sols, nonzero_ids_unknowns[j], :)
            md_sols_part = md[1:n_unknowns][nonzero_ids_unknowns[j]]
            evaluated_mons[j, :, i] = vec(prod(sols_part.^md_sols_part, dims=1)).*params_eval
        end
    end
    return evaluated_mons
end