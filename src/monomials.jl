export Monomial,
    DenseMonomialVector,
    multiexponents,
    iterate,
    extend!,
    evaluate,
    to_expressions,
    to_classes,
    only_param_dep,
    n_only_param_dep

# TODO: remove?
struct Monomial{T<:Integer,S<:Integer}
    unkn_mexp::SparseVector{T,S}
    param_mexp::SparseVector{T,S}
    unknowns::Vector{Variable}
    parameters::Vector{Variable}

    function Monomial{T,S}(
        unkn_mexp,
        param_mexp,
        unknowns,
        parameters
    ) where {T,S<:Integer}
        # TODO: check if lengths are equal, if mds contains numbers >= 0
        return new(unkn_mexp, param_mexp, unknowns, parameters)
    end
end

Monomial(
    md::SparseVector{T,S},
    vars::Vector{Variable}
) where {T,S<:Integer} = Monomial{T,S}(md, vars)

function Monomial{T,S}(var::Variable, vars::Vector{Variable}) where {T,S<:Integer}
    md = spzeros(T, S, length(vars))
    md[findfirst(x->x==var, vars)] = 1
    return Monomial{T,S}(md, vars)
end

function Monomial{T,S}(
    mon::Expression,
    vars::Vector{Variable}
) where {T,S<:Integer}
    es, cs = exponents_coefficients(mon, vars)
    if length(cs) > 1 || !isone(cs[1])
        throw(ArgumentError("Input expression is not a monomial"))
    else
        return Monomial{T,S}(sparse(es[:,1]), vars)
    end
end

prodpow(v::AbstractVector, e::AbstractSparseVector) = prod(v[e.nzind].^e.nzval)
Base.isone(mon::Monomial) = iszero(mon.md)
Base.convert(
    ::Type{Expression},
    mon::Monomial
) = prodpow(mon.unknowns, mon.unkn_mexp)*prodpow(mon.parameters, mon.param_mexp)
Base.:(==)(m₁::Monomial, m₂::Monomial) = Expression(m₁) == Expression(m₂)
Base.show(io::IO, mon::Monomial) = show(io, Expression(mon))


"""
    multiexponents(; degree::T, nvars::S) where {T<:Integer,S<:Integer}

Generates vector of multiexponents of total degree `degree` each represented by `SparseVector{T,S}` of length `nvars`.
"""
function multiexponents(; degree::T, nvars::S) where {T<:Integer,S<:Integer}
    mexps = [spzeros(T, S, nvars) for _ in 1:num_mons(nvars, degree)]
    k = 1
    for n in 1:nvars
        for part::Vector{T} in partitions(degree, n)
            for vals in multiset_permutations(part, n)
                for inds in combinations(S.(1:nvars), n)
                    mexps[k][inds] = vals
                    k += 1
                end
            end
        end
    end
    return mexps
end

abstract type AbstractMonomialVector end

# Structure representing vector of monomials
# Advantage: monomial evalution can be made more efficient
mutable struct DenseMonomialVector{T<:Integer,S<:Integer} <: AbstractMonomialVector
    degree::T
    unknowns_mexps::Dict{T, Vector{SparseVector{T,S}}}  # keys are degrees
    parameters_mexps::Dict{T, Vector{SparseVector{T,S}}}  # keys are degrees
    unknowns::Vector{Variable}
    parameters::Vector{Variable}
end

function DenseMonomialVector{T,S}(;
    unknowns::Vector{Variable},
    parameters::Vector{Variable}
) where {T,S<:Integer}

    degree = zero(T)
    unknowns_mexps = Dict(0 => [spzeros(T, S, length(unknowns))])
    parameters_mexps = Dict(0 => [spzeros(T, S, length(parameters))])
    return DenseMonomialVector{T,S}(
        degree,
        unknowns_mexps,
        parameters_mexps,
        unknowns,
        parameters
    )
end

unknowns(mons::DenseMonomialVector) = mons.unknowns
HC.parameters(mons::DenseMonomialVector) = mons.parameters
variables(mons::DenseMonomialVector) = vcat(unknowns(mons), parameters(mons))

nunknowns(mons::DenseMonomialVector{T,S}) where {T<:Integer,S<:Integer} = S(length(mons.unknowns))
HC.nparameters(mons::DenseMonomialVector{T,S}) where {T<:Integer,S<:Integer} = S(length(mons.parameters))
nvariables(mons::DenseMonomialVector{T,S}) where {T<:Integer,S<:Integer} = nunknowns(mons)+nparameters(mons)

Base.length(mons::DenseMonomialVector) = num_mons_upto(nvariables(mons), mons.degree)
nparam_only(mons::DenseMonomialVector) = num_mons_upto(nparameters(mons), mons.degree)
nunkn_dep(mons::DenseMonomialVector) = length(mons) - nparam_only(mons)

to_expression(
    unknowns::Vector{Variable},
    unkn_mexp::AbstractSparseVector,
    parameters::Vector{Variable},
    param_mexp::AbstractSparseVector
) = prodpow(unknowns, unkn_mexp)*prodpow(parameters, param_mexp)

function to_expressions(mons::DenseMonomialVector)
    return [to_expression(
                mons.unknowns,
                unkn_mexp,
                mons.parameters,
                param_mexp
            ) for (unkn_mexp, param_mexp) in mons]
end

function Base.show(io::IO, mons::DenseMonomialVector)
    println(io, "$(length(mons))-element $(typeof(mons))")
    print(io, "[", join(to_expressions(mons), ", "), "]")
end

"""
    extend!(mons::DenseMonomialVector{T,S}; degree::Integer) where {T<:Integer,S<:Integer}

Extends monomial vector `mons` by new multiexponents of total degree upto `degree`.
"""
function extend!(mons::DenseMonomialVector{T,S}; degree::Integer) where {T<:Integer,S<:Integer}
    degree = convert(T, degree)
    if degree > mons.degree
        for d::T in mons.degree+1:degree
            mons.unknowns_mexps[d] = multiexponents(degree=d, nvars=nunknowns(mons))
            mons.parameters_mexps[d] = multiexponents(degree=d, nvars=nparameters(mons))
        end
        mons.degree = degree
    end
    return nothing
end

struct DenseMonVecState
    udeg::Int
    uid::Int
    pdeg::Int
    pid::Int
end

function pick_at_state(
    mons::DenseMonomialVector,
    state::DenseMonVecState
)
    @unpack udeg, uid, pdeg, pid = state
    uexps, pexps = mons.unknowns_mexps, mons.parameters_mexps
    return uexps[udeg][uid], pexps[pdeg][pid]
end

function next_state(
    mons::DenseMonomialVector,
    state::DenseMonVecState
)
    @unpack udeg, uid, pdeg, pid = state
    uexps, pexps = mons.unknowns_mexps, mons.parameters_mexps
    if uid < length(uexps[udeg])
        return DenseMonVecState(udeg, uid+1, pdeg, pid)
    elseif pid < length(pexps[pdeg])
        return DenseMonVecState(udeg, 1, pdeg, pid+1)
    elseif pdeg > 0
        return DenseMonVecState(udeg, 1, pdeg-1, 1)
    elseif udeg > 0
        return DenseMonVecState(udeg-1, 1, mons.degree-udeg+1, 1)
    else
        return nothing
    end
end

"""
    iterate(mons::DenseMonomialVector, state::DenseMonVecState)

Returns the next item (multiexponent) together with the next state. 
The state is defined to be a tuple of (udeg, uid, pdeg, pid), where 
`udeg` and `pdeg` are the keys for `mons.unknowns_mexps` and `mons.parameters_mexps`,
respectively, and `uid` and `pid` are the indices in the corresponding values of 
these dictionaries.
"""
function Base.iterate(
    mons::DenseMonomialVector,
    state::DenseMonVecState
)
    new_state = next_state(mons, state)
    isnothing(new_state) && return nothing
    return pick_at_state(mons, new_state), new_state
end

"""
    iterate(mons::DenseMonomialVector)

Returns the first item (multiexponent) together with the state.
"""
function Base.iterate(mons::DenseMonomialVector)
    state = DenseMonVecState(mons.degree, 1, 0, 1)
    return pick_at_state(mons, state), state
end

# Structure representing sparse monomial vector
# Advantages:
#   1) monomial evalution can be made more efficient
#   2) we know how many param_only monomials are there (needed for picking samples for evaluation)
struct SparseMonomialVector{T<:Integer,S<:Integer} <: AbstractMonomialVector
    unkn_dep_mexps::NTuple{2, Vector{SparseVector{T,S}}}
    param_only_mexps::Vector{SparseVector{T,S}}
    unknowns::Vector{Variable}
    parameters::Vector{Variable}
end

SparseMonomialVector{Tv,Ti}(
    unknowns::Vector{Variable},
    parameters::Vector{Variable}
) where {Tv<:Integer,Ti<:Integer} = SparseMonomialVector{Tv,Ti}(([], []), [], unknowns, parameters)

Base.length(mons::SparseMonomialVector) = length(mons.unkn_dep_mexps[1]) + length(mons.param_only_mexps)
is_param_only(mons::SparseMonomialVector) = length(mons.unkn_dep_mexps[1]) == 0
nunkn_dep(mons::SparseMonomialVector) = length(mons.unkn_dep_mexps[1])
nparam_only(mons::SparseMonomialVector) = length(mons.param_only_mexps)

function Base.push!(
    mons::SparseMonomialVector{T,S},
    (uexp, pexp)::NTuple{2, SparseVector{T,S}}
) where {T<:Integer,S<:Integer}

    if iszero(uexp)
        push!(mons.param_only_mexps, pexp)
    else
        push!(mons.unkn_dep_mexps[1], uexp)
        push!(mons.unkn_dep_mexps[2], pexp)
    end
end

function to_expressions(mons::SparseMonomialVector)
    return vcat(
        [prodpow(mons.unknowns, unkn_md)*prodpow(mons.parameters, param_md) for (unkn_md, param_md) in mons.mixed_mds],
        [prodpow(mons.parameters, md) for md in mons.param_only_mexps]
    )
end

# TODO
function Base.show(io::IO, mons::SparseMonomialVector)
    # println(io, "$(length(mons.mds))-element $(typeof(mons))")
    print(io, "[", join(to_expressions(mons), ", "), "]")
end

function Base.gcd(mons::SparseMonomialVector)
    return Monomial(
        min.(mons.unkn_dep_mexps[1]...),
        min.(mons.unkn_dep_mexps[2]..., mons.param_only_mexps...),
        mons.unknowns,
        mons.parameters
    )
end

# Structure for evaluated monomials
# Advatage: requires less memory (no need to duplicate values for param_only monomials)
struct EvaluatedMonomials{T <: AbstractMonomialVector}
    unkn_dep::Array{ComplexF64, 3}
    param_only::Array{ComplexF64, 2}
    mons::T
end

# TODO: test timings w/ and w/o selectdim
prodpow(
    v::AbstractArray,
    e::AbstractSparseVector;
) = dropdims(prod(selectdim(v,1,e.nzind).^e.nzval; dims=1); dims=1)

function HC.evaluate(
    mexps_dict::Dict{T, Vector{SparseVector{T,S}}},
    samples::Array{P,N}
) where {T, S, P, N}

    evals = Dict{T, Array{P,N}}()
    for (d, mexps) in mexps_dict
        eval = zeros(P, length(mexps), size(samples)[2:end]...)
        for (i, mexp) in enumerate(mexps)
            eval[i, repeat([:], N-1)...] = prodpow(samples, mexp)
        end
        evals[d] = eval
    end
    return evals
end

function HC.evaluate(
    mons::DenseMonomialVector,
    samples::Samples
)
    unkn_evals = evaluate(mons.unknowns_mexps, samples.solutions)
    param_evals = evaluate(mons.parameters_mexps, samples.parameters)
    unkn_dep = zeros(ComplexF64, nunkn_dep(mons), nsolutions(samples), ninstances(samples))
    param_only = zeros(ComplexF64, nparam_only(mons), ninstances(samples))
    k = 0
    for udeg in mons.degree:-1:1
        unkn_eval = unkn_evals[udeg]
        for pdeg in (mons.degree-udeg):-1:0
            param_eval = param_evals[pdeg]
            for pid in axes(param_eval, 1)
                unkn_dep[(k+1):(k+size(unkn_eval, 1)), :, :] = unkn_eval.*reshape(param_eval[pid, :], 1, 1, :)
                k += size(unkn_eval, 1)
            end
        end
    end
    param_only = vcat([param_evals[pdeg] for pdeg in mons.degree:-1:0]...)
    return EvaluatedMonomials(unkn_dep, param_only, mons)
end

function HC.evaluate(
    mons::SparseMonomialVector,
    samples::Samples
)
    param_only = zeros(ComplexF64, nparam_only(mons), ninstances(samples))
    for (i, mexp) in enumerate(mons.param_only_mexps)
        param_only[i, :] = prodpow(samples.parameters, mexp)
    end
    unkn_dep = zeros(ComplexF64, nmixed(mons), nsolutions(samples), ninstances(samples))
    for (i, (uexp, pexp)) in enumerate(zip(mons.unkn_dep_mexps...))
        unkn_dep[i, :, :] = prodpow(samples.solutions, uexp).*prodpow(samples.parameters, pexp)'
    end
    return EvaluatedMonomials(unkn_dep, param_only, mons)
end

HC.evaluate(
    mons::AbstractMonomialVector,
    samples::Dict{Vector{Int}, Samples}
) = Dict(zip(keys(samples), [evaluate(mons, s) for (_, s) in samples]))
