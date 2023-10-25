export ScalingGroup,
    scaling_symmetries

using HomotopyContinuation: exponents_coefficients
using AbstractAlgebra: matrix, snf_with_transform, ZZ, hnf, GF, lift

Grading = Vector{Tuple{Int, Matrix{Int}}}  # TODO: make it struct?

struct ScalingGroup
    grading::Grading
    vars::Vector{Variable}
    action::Vector{Tuple{Int, Dict{Variable, Expression}}}
end

ScalingGroup() = ScalingGroup([], [], [])

# TODO: replace finite λ vars to multiplications by roots of unity?
function ScalingGroup(grading::Grading, vars::Vector{Variable})
    action = Vector{Tuple{Int, Dict{Variable, Expression}}}([])
    k = 0
    for (sᵢ, Uᵢ) in grading
        n_scalings = size(Uᵢ, 1)
        @var λ[(k+1):(n_scalings+k)]
        for j in 1:n_scalings
            nonzero_ids = findall(!iszero, Uᵢ[j, :])
            vals = (λ[j].^Uᵢ[j, :][nonzero_ids]).*vars[nonzero_ids]
            push!(action, (sᵢ, Dict(zip(vars[nonzero_ids], vals))))
        end
        k += n_scalings
    end
    return ScalingGroup(grading, vars, action)
end

Base.copy(s::ScalingGroup) = ScalingGroup(s.grading, s.vars, s.action)

function Base.show(io::IO, scalings::ScalingGroup)
    n_infinite, n_finite = 0, 0
    for (sᵢ, Uᵢ) in scalings.grading
        if sᵢ == 0
            n_infinite = size(Uᵢ, 1)
        else
            n_finite += size(Uᵢ, 1)
        end
    end
    print(io, "ScalingGroup with $(n_infinite+n_finite) scaling")
    n_infinite + n_finite == 1 ? print(io, "\n") : print(io, "s\n")
    n_infinite + n_finite == 0 && return
    if n_infinite != 0
        print(io, " $(n_infinite) infinite scaling")
        n_infinite == 1 ? print(io, ":") : print(io, "s:")
        for (sᵢ, scaling) in scalings.action
            if sᵢ == 0
                print(io, "\n  ")
                for (j, var) in enumerate(scalings.vars)
                    if !isnothing(get(scaling, var, nothing))
                        print(io, var, " ↦ ", scaling[var])
                        j < length(scalings.vars) && print(io, ", ")
                    end
                end
            end
        end
    end
    if n_finite != 0
        println(io, " finite scalings:")
        for (i, (sᵢ, Uᵢ)) in enumerate(scalings.grading)
            if sᵢ == 0 continue end
            println(io, "  $(size(Uᵢ, 1)) of order $(sᵢ)")
        end
    end
    # print(io, " vars: ", join(scalings.vars, ", "))
end

function _mat2col_diffs(M::AbstractMatrix{T}) where {T<:Number}
    M = M - M[:,1]*ones(T, 1, size(M,2))
    return M[:,2:end]
end

function _snf_scaling_symmetries(F::System)::Tuple{Vector{Int}, Vector{Matrix{Int}}}
    vars = vcat(F.variables, F.parameters)
    Es = [exponents_coefficients(f, vars)[1] for f in F.expressions]
    K = hcat([_mat2col_diffs(E) for E in Es]...)
    if size(K, 1) > size(K, 2)
        K = [K zeros(eltype(K), size(K, 1), size(K,1)-size(K,2))]
    end
    K = matrix(ZZ, K)

    S, U, _ = snf_with_transform(K)
    U, S = Matrix(U), Int.(diag(Matrix(S)))

    s = reverse(filter(el->el!=1, unique(S)))
    if length(s) == 0
        return [], []
    end
    idxs = [findall(x->x==el, S) for el in s]
    Us = [U[idxs[i], :] for i in eachindex(idxs)]
    return s, Us # TODO: what if max(Us[i]) > MAX_INT64?
end

function _hnf_reduce!(grading::Grading)
    for (i, (sᵢ, Uᵢ)) in enumerate(grading)
        if sᵢ == 0
            grading[i] = (sᵢ, Matrix(hnf(matrix(ZZ, Uᵢ))))
        else
            grading[i] = (sᵢ, lift.(Matrix(hnf(matrix(GF(sᵢ), Uᵢ)))))
        end
    end
end

"""
    scaling_symmetries(F::System; in_hnf::Bool=true)

Given a polynomial system `F` returns the group of scaling symmetries 
of the polynomial system `F`.

```julia-repl
julia> @var x[1:2] p[1:2];

julia> F = System([x[1]^2 - x[2]^2 - p[1], 2*x[1]*x[2] - p[2]]; variables=x, parameters=p);

julia> scalings = scaling_symmetries(F)
ScalingGroup with 2 scalings
 infinite scalings: 1
 finite scalings:
  1 of order 2
 vars: x₁, x₂, p₁, p₂
```
"""
function scaling_symmetries(F::System; in_hnf::Bool=true)
    s, U = _snf_scaling_symmetries(F)
    if length(s) == 0
        return ScalingGroup()
    end
    grading = collect(zip(s, U))
    if in_hnf _hnf_reduce!(grading) end
    vars = vcat(F.variables, F.parameters)
    return ScalingGroup(grading, vars)
end

scaling_symmetries(
    F::SampledSystem;
    in_hnf::Bool=true
) = scaling_symmetries(F.system; in_hnf=in_hnf)

# TODO: extend to remove rows dependent on other blocks
function reduce(grading::Grading)
    grading_cp = copy(grading)
    _hnf_reduce!(grading_cp)
    filtered_grading = Grading([])
    for (i, (sᵢ, Uᵢ)) in enumerate(grading_cp)
        Uᵢ = remove_zero_rows(Uᵢ)
        if size(Uᵢ, 1) != 0
            push!(filtered_grading, (sᵢ, Uᵢ))
        end
    end
    return filtered_grading
end

function restrict_scalings(scalings::ScalingGroup, var_ids::Vector{Int})
    restr_grading = copy(scalings.grading)
    for (i, (sᵢ, Uᵢ)) in enumerate(scalings.grading)
        restr_grading[i] = (sᵢ, Uᵢ[:, idx])
    end
    return ScalingGroup(reduce(restr_grading), scalings.vars[var_ids])
end

# TODO: what if vars is not a subset of scalings.vars?
function restrict_scalings(scalings::ScalingGroup, vars::Vector{Variable})
    ids = [findfirst(v->v==var, scalings.vars) for var in vars]
    if nothing in ids 
        throw(ArgumentError("vars contains variables not present in scalings.vars"))
    end
    return restrict_scalings(scalings, ids)
end

scaling_symmetries(
    F::System,
    vars::Vector{Variable}
) = restrict_scalings(scaling_symmetries(F), vars)
scaling_symmetries(
    F::SampledSystem,
    vars::Vector{Variable}
) = scaling_symmetries(F.system, vars)

function HC.degree(md::Vector{<:Integer}, grading::Grading)
    return vcat([mod(Uᵢ*md, sᵢ) for (sᵢ, Uᵢ) in grading]...)
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