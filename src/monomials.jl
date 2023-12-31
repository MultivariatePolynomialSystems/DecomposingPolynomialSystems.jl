export multiexponents,
    AbstractMonomialVector,
    DenseMonomialVector,
    SparseMonomialVector,
    iterate,
    extend!,
    evaluate,
    to_expressions,
    to_classes,
    nparam_only,
    nunkn_dep

function multiexponents(; degree::Tv, nvars::Ti) where {Tv<:Integer,Ti<:Integer}
    mexps = [spzeros(Tv, Ti, nvars) for _ in 1:num_mons(nvars, degree)]
    k = 1
    for n in 1:nvars
        for part::Vector{Tv} in partitions(degree, n)
            for vals in multiset_permutations(part, n)
                for inds in combinations(Ti.(1:nvars), n)
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
mutable struct DenseMonomialVector{Tv<:Integer,Ti<:Integer} <: AbstractMonomialVector
    degree::Tv
    n::Int  # length
    unknowns_mexps::Dict{Tv, Vector{SparseVector{Tv,Ti}}}  # keys are degrees
    parameters_mexps::Dict{Tv, Vector{SparseVector{Tv,Ti}}}  # keys are degrees
    unknowns::Vector{Variable}
    parameters::Vector{Variable}
end

function DenseMonomialVector{Tv,Ti}(;
    unknowns::Vector{Variable},
    parameters::Vector{Variable}
) where {Tv<:Integer,Ti<:Integer}

    degree = zero(Tv)
    unknowns_mexps = Dict(0 => [spzeros(Tv, Ti, length(unknowns))])
    parameters_mexps = Dict(0 => [spzeros(Tv, Ti, length(parameters))])
    n = 1
    return DenseMonomialVector{Tv,Ti}(
        degree,
        n,
        unknowns_mexps,
        parameters_mexps,
        unknowns,
        parameters
    )
end

unknowns(mons::DenseMonomialVector) = mons.unknowns
HC.parameters(mons::DenseMonomialVector) = mons.parameters
variables(mons::DenseMonomialVector) = vcat(unknowns(mons), parameters(mons))

nunknowns(mons::DenseMonomialVector{Tv,Ti}) where {Tv<:Integer,Ti<:Integer} = length(mons.unknowns)
HC.nparameters(mons::DenseMonomialVector{Tv,Ti}) where {Tv<:Integer,Ti<:Integer} = length(mons.parameters)
nvariables(mons::DenseMonomialVector{Tv,Ti}) where {Tv<:Integer,Ti<:Integer} = nunknowns(mons)+nparameters(mons)

unknowns_deg(mons::DenseMonomialVector) = max(keys(mons.unknowns_mexps)...)
parameters_deg(mons::DenseMonomialVector) = max(keys(mons.parameters_mexps)...)

Base.length(mons::DenseMonomialVector) = mons.n
nparam_only(mons::DenseMonomialVector) = num_mons_upto(nparameters(mons), parameters_deg(mons))
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

function extend!(
    mons::DenseMonomialVector{Tv,Ti};
    degree::Integer,
    extend_params::Bool=true
) where {Tv<:Integer,Ti<:Integer}
    degree = convert(Tv, degree)
    uexps, pexps = mons.unknowns_mexps, mons.parameters_mexps
    for d::Tv in mons.degree+1:degree
        uexps[d] = multiexponents(degree=d, nvars=Ti(nunknowns(mons)))
    end
    if extend_params
        for d::Tv in min(keys(pexps)...)+1:degree
            pexps[d] = multiexponents(degree=d, nvars=Ti(nparameters(mons)))
        end
    end
    mons.degree = degree
    mons.n = 0
    for udeg in 0:unknowns_deg(mons)
        for pdeg in 0:parameters_deg(mons)
            if udeg+pdeg ≤ mons.degree
                mons.n += num_mons(nunknowns(mons), udeg)*num_mons(nparameters(mons), pdeg)
            else
                break
            end
        end
    end
    return mons
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
        return DenseMonVecState(udeg-1, 1, min(mons.degree-udeg+1, parameters_deg(mons)), 1)
    else
        return nothing
    end
end

function Base.iterate(
    mons::DenseMonomialVector,
    state::DenseMonVecState
)
    new_state = next_state(mons, state)
    isnothing(new_state) && return nothing
    return pick_at_state(mons, new_state), new_state
end

function Base.iterate(mons::DenseMonomialVector)
    state = DenseMonVecState(mons.degree, 1, 0, 1)
    return pick_at_state(mons, state), state
end

# Structure representing sparse monomial vector
# Advantages:
#   1) monomial evalution can be made more efficient
#   2) we know how many param_only monomials are there (needed for picking samples for evaluation)
struct SparseMonomialVector{Tv<:Integer,Ti<:Integer} <: AbstractMonomialVector
    unkn_dep_mexps::NTuple{2, Vector{SparseVector{Tv,Ti}}}
    param_only_mexps::Vector{SparseVector{Tv,Ti}}
    unknowns::Vector{Variable}
    parameters::Vector{Variable}
end

SparseMonomialVector{Tv,Ti}(
    unknowns::Vector{Variable},
    parameters::Vector{Variable}
) where {Tv<:Integer,Ti<:Integer} = SparseMonomialVector{Tv,Ti}(([], []), [], unknowns, parameters)

struct SparseMonVecState
    unkn_dep::Bool
    id::Int
end

function pick_at_state(mons::SparseMonomialVector{Tv,Ti}, state::SparseMonVecState) where {Tv<:Integer,Ti<:Integer}
    state.unkn_dep && return vcat(mons.unkn_dep_mexps[1][state.id], mons.unkn_dep_mexps[2][state.id])
    return vcat(spzeros(Tv,Ti,nunknowns(mons)), mons.param_only_mexps[state.id])
end

function Base.iterate(mons::SparseMonomialVector)
    length(mons) == 0 && return nothing
    unkn_dep = mons.unkn_dep_mexps
    if length(unkn_dep[1]) != 0
        state = SparseMonVecState(true, 1)
        return pick_at_state(mons, state), state
    else
        state = SparseMonVecState(false, 1)
        return pick_at_state(mons, state), state
    end
end

function next_state(
    mons::SparseMonomialVector,
    state::SparseMonVecState
)
    @unpack unkn_dep, id = state
    if unkn_dep
        if id < length(mons.unkn_dep_mexps[1])
            return SparseMonVecState(true, id+1)
        else
            if length(mons.param_only_mexps) != 0
                return SparseMonVecState(false, 1)
            else
                return nothing
            end
        end
    else
        if id < length(mons.param_only_mexps)
            return SparseMonVecState(true, id+1)
        else
            return nothing
        end
    end
end

function Base.iterate(
    mons::SparseMonomialVector,
    state::SparseMonVecState
)
    new_state = next_state(mons, state)
    isnothing(new_state) && return nothing
    return pick_at_state(mons, new_state), new_state
end

function Base.findfirst(mexp::SparseVector{Tv,Ti}, mons::SparseMonomialVector{Tv,Ti}) where {Tv<:Integer,Ti<:Integer}
    for (i, mon) in enumerate(mons)
        mexp == mon && return i
    end
    return nothing
end

unknowns(mons::SparseMonomialVector) = mons.unknowns
HC.parameters(mons::SparseMonomialVector) = mons.parameters
variables(mons::SparseMonomialVector) = vcat(unknowns(mons), parameters(mons))
nunknowns(mons::SparseMonomialVector) = length(mons.unknowns)
HC.nparameters(mons::SparseMonomialVector) = length(mons.parameters)
Base.length(mons::SparseMonomialVector) = length(mons.unkn_dep_mexps[1]) + length(mons.param_only_mexps)
is_param_only(mons::SparseMonomialVector) = length(mons.unkn_dep_mexps[1]) == 0
nunkn_dep(mons::SparseMonomialVector) = length(mons.unkn_dep_mexps[1])
nparam_only(mons::SparseMonomialVector) = length(mons.param_only_mexps)

function Base.vcat(monVs::SparseMonomialVector{Tv,Ti}...) where {Tv<:Integer,Ti<:Integer}
    return SparseMonomialVector{Tv,Ti}(
        (vcat([mons.unkn_dep_mexps[1] for mons in monVs]...), vcat([mons.unkn_dep_mexps[2] for mons in monVs]...)),
        vcat([mons.param_only_mexps for mons in monVs]...),
        monVs[1].unknowns,
        monVs[1].parameters
    )
end

function Base.push!(
    mons::SparseMonomialVector{Tv,Ti},
    (uexp, pexp)::NTuple{2, SparseVector{Tv,Ti}}
) where {Tv<:Integer,Ti<:Integer}

    if iszero(uexp)
        push!(mons.param_only_mexps, pexp)
    else
        push!(mons.unkn_dep_mexps[1], uexp)
        push!(mons.unkn_dep_mexps[2], pexp)
    end
end

function to_expressions(mons::SparseMonomialVector)
    return vcat(
        [to_expression(
            mons.unknowns,
            unkn_mexp,
            mons.parameters,
            param_mexp
        ) for (unkn_mexp, param_mexp) in zip(mons.unkn_dep_mexps...)],
        [prodpow(mons.parameters, mexp) for mexp in mons.param_only_mexps]
    )
end

# TODO
function Base.show(io::IO, mons::SparseMonomialVector)
    # println(io, "$(length(mons.mds))-element $(typeof(mons))")
    print(io, "[", join(to_expressions(mons), ", "), "]")
end

function Base.gcd(mons::SparseMonomialVector{Tv,Ti}) where {Tv<:Integer,Ti<:Integer}
    if nparam_only(mons) == 0
        return vcat(
            min.(mons.unkn_dep_mexps[1]...),
            min.(mons.unkn_dep_mexps[2]...)
        )
    end
    return vcat(
        spzeros(Tv, Ti, nunknowns(mons)),
        min.(mons.unkn_dep_mexps[2]..., mons.param_only_mexps...)
    )
end

struct SparseMatrixMonomialVector{Tv<:Integer,Ti<:Integer} <: AbstractMonomialVector
    unkn_dep_mexps::NTuple{2, SparseMatrixCSC{Tv,Ti}}
    param_only_mexps::SparseMatrixCSC{Tv,Ti}
    unknowns::Vector{Variable}
    parameters::Vector{Variable}
end

function SparseMatrixMonomialVector{Tv,Ti}(
    mons::SparseMonomialVector{Tv,Ti}
) where {Tv<:Integer,Ti<:Integer}

    unkn_dep_1 = length(mons.unkn_dep_mexps[1]) == 0 ? spzeros(nunknowns(mons), 0) : Base.reduce(hcat, mons.unkn_dep_mexps[1])
    unkn_dep_2 = length(mons.unkn_dep_mexps[2]) == 0 ? spzeros(nparameters(mons), 0) : Base.reduce(hcat, mons.unkn_dep_mexps[2])
    param_only = length(mons.param_only_mexps) == 0 ? spzeros(nparameters(mons), 0) : Base.reduce(hcat, mons.param_only_mexps)
    return SparseMatrixMonomialVector{Tv,Ti}(
        (unkn_dep_1, unkn_dep_2),
        param_only,
        mons.unknowns,
        mons.parameters
    )
end

nunknowns(mons::SparseMatrixMonomialVector) = length(mons.unknowns)
Base.length(mons::SparseMatrixMonomialVector) = size(mons.unkn_dep_mexps[1], 2) + size(mons.param_only_mexps, 2)
is_param_only(mons::SparseMatrixMonomialVector) = size(mons.unkn_dep_mexps[1], 2) == 0
nunkn_dep(mons::SparseMatrixMonomialVector) = size(mons.unkn_dep_mexps[1], 2)
nparam_only(mons::SparseMatrixMonomialVector) = size(mons.param_only_mexps, 2)

function Base.vcat(monVs::SparseMatrixMonomialVector{Tv,Ti}...) where {Tv<:Integer,Ti<:Integer}
    return SparseMatrixMonomialVector{Tv,Ti}(
        (Base.reduce(hcat, [mons.unkn_dep_mexps[1] for mons in monVs]), Base.reduce(hcat, [mons.unkn_dep_mexps[2] for mons in monVs])),
        Base.reduce(hcat, [mons.param_only_mexps for mons in monVs]),
        monVs[1].unknowns,
        monVs[1].parameters
    )
end

function to_expressions(mons::SparseMatrixMonomialVector)
    return vcat(
        [to_expression(
            mons.unknowns,
            mons.unkn_dep_mexps[1][:, i],
            mons.parameters,
            mons.unkn_dep_mexps[2][:, i]
        ) for i in axes(mons.unkn_dep_mexps[1], 2)],
        [prodpow(mons.parameters, mons.param_only_mexps[:, i]) for i in axes(mons.param_only_mexps, 2)]
    )
end

function Base.gcd(mons::SparseMatrixMonomialVector{Tv,Ti}) where {Tv<:Integer,Ti<:Integer}
    if nparam_only(mons) == 0
        return vcat(
            vec(minimum(mons.unkn_dep_mexps[1]; dims=2)),
            vec(minimum(mons.unkn_dep_mexps[2]; dims=2))
        )
    elseif nunkn_dep(mons) == 0
        return vcat(
            spzeros(Tv, Ti, nunknowns(mons)),
            vec(minimum(mons.param_only_mexps; dims=2))
        )
    else
        return vcat(
            spzeros(Tv, Ti, nunknowns(mons)),
            vec(minimum(hcat(minimum(mons.unkn_dep_mexps[2]; dims=2), minimum(mons.param_only_mexps; dims=2)); dims=2))
        )
    end
end

# Structure for evaluated monomials
# Advatage: requires less memory (no need to duplicate values for param_only monomials)
struct EvaluatedMonomials
    unkn_dep::Array{ComplexF64, 3}
    param_only::Array{ComplexF64, 2}
end

nunkn_dep(eval_mons::EvaluatedMonomials) = size(eval_mons.unkn_dep, 1)
nparam_only(eval_mons::EvaluatedMonomials) = size(eval_mons.param_only, 1)
nmonomials(eval_mons::EvaluatedMonomials) = nunkn_dep(eval_mons) + nparam_only(eval_mons)

# TODO: test timings w/ and w/o selectdim
prodpow(
    v::AbstractArray,
    e::AbstractSparseVector;
) = dropdims(prod(selectdim(v,1,e.nzind).^e.nzval; dims=1); dims=1)

function HC.evaluate(
    mexps_dict::Dict{Tv, Vector{SparseVector{Tv,Ti}}},
    samples::Array{P,N}
) where {Tv<:Integer, Ti<:Integer, P, N}

    evals = Dict{Tv, Array{P,N}}()
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
    for udeg in unknowns_deg(mons):-1:1
        unkn_eval = unkn_evals[udeg]
        for pdeg in min(mons.degree-udeg, parameters_deg(mons)):-1:0
            param_eval = param_evals[pdeg]
            for pid in axes(param_eval, 1)
                unkn_dep[(k+1):(k+size(unkn_eval, 1)), :, :] = unkn_eval.*reshape(param_eval[pid, :], 1, 1, :)
                k += size(unkn_eval, 1)
            end
        end
    end
    param_only = vcat([param_evals[pdeg] for pdeg in parameters_deg(mons):-1:0]...)  # TODO: replace vcat?
    return EvaluatedMonomials(unkn_dep, param_only)
end

function HC.evaluate(
    mons::SparseMonomialVector,
    samples::Samples
)
    param_only = zeros(ComplexF64, nparam_only(mons), ninstances(samples))
    for (i, mexp) in enumerate(mons.param_only_mexps)
        param_only[i, :] = prodpow(samples.parameters, mexp)
    end
    unkn_dep = zeros(ComplexF64, nunkn_dep(mons), nsolutions(samples), ninstances(samples))
    for (i, (uexp, pexp)) in enumerate(zip(mons.unkn_dep_mexps...))
        unkn_dep[i, :, :] = prodpow(samples.solutions, uexp).*transpose(prodpow(samples.parameters, pexp))
    end
    return EvaluatedMonomials(unkn_dep, param_only)
end

function HC.evaluate(
    mons::SparseMatrixMonomialVector,
    samples::Samples
)
    param_only = zeros(ComplexF64, nparam_only(mons), ninstances(samples))
    for i in axes(mons.param_only_mexps, 2)
        param_only[i, :] = prodpow(samples.parameters, mons.param_only_mexps[:, i])
    end
    unkn_dep = zeros(ComplexF64, nunkn_dep(mons), nsolutions(samples), ninstances(samples))
    for i in axes(mons.unkn_dep_mexps[1], 2)
        unkn_dep[i, :, :] = prodpow(samples.solutions, mons.unkn_dep_mexps[1][:, i]).*transpose(prodpow(samples.parameters, mons.unkn_dep_mexps[2][:, i]))
    end
    return EvaluatedMonomials(unkn_dep, param_only)
end

HC.evaluate(
    mons::AbstractMonomialVector,
    samples::Dict{T, Samples}
) where {T} = Dict(zip(values(samples), [evaluate(mons, s) for s in values(samples)]))
