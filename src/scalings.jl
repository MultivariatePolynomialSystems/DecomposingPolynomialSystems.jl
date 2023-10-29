export Grading,
    ScalingGroup,
    scaling_symmetries

using AbstractAlgebra: matrix, snf_with_transform, ZZ, hnf, GF, lift
using LinearAlgebra: diag

mutable struct Grading
    free_part::Union{Nothing, Matrix{Int}}
    mod_part::Vector{Tuple{Int, Matrix{Int}}}
end

Grading() = Grading(nothing, [])

function Grading(s::Vector{Int}, U::Vector{Matrix{Int}})
    grading = Grading()
    for (sᵢ, Uᵢ) in zip(s, U)
        if sᵢ == 0
            grading.free_part = Uᵢ
        else
            push!(grading.mod_part, (sᵢ, Uᵢ))
        end
    end
    return grading
end

Base.copy(grading::Grading) = Grading(grading.free_part, grading.mod_part)

function _structure(grading::Grading)
    str = ""
    U₀ = grading.free_part
    if !isnothing(U₀)
        n_free = size(U₀, 1)
        free_str = n_free == 1 ? "ℤ" : "ℤ$(superscript(n_free))"
        str = str * free_str
    end
    for (i, (sᵢ, Uᵢ)) in enumerate(grading.mod_part)
        prod_str = isnothing(U₀) && i == 1 ? "" : " × "
        n_sᵢ = size(Uᵢ, 1)
        sᵢ_str = n_sᵢ == 1 ? "ℤ$(subscript(sᵢ))" : "ℤ$(subscript(sᵢ))$(superscript(n_sᵢ))"
        str = str * prod_str * sᵢ_str
    end
    return str
end

function Base.show(io::IO, grading::Grading)

end

SparseAction = Vector{Tuple{Variable, Expression}}

struct ScalingGroup
    grading::Grading
    structure::String
    vars::Vector{Variable}
    action::Tuple{Vector{SparseAction}, Vector{Tuple{Int, Vector{SparseAction}}}}
end

ScalingGroup(vars::Vector{Variable}) = ScalingGroup(Grading(), "N/A", vars, ([], []))

function ScalingGroup(grading::Grading, vars::Vector{Variable})
    free_action = Vector{SparseAction}([])
    U₀ = grading.free_part
    if !isnothing(U₀)
        if size(U₀, 1) == 1
            @var λ
            λ = [λ]
        else
            @var λ[1:size(U₀, 1)]
        end
        for j in axes(U₀, 1)
            nonzero_ids = findall(!iszero, U₀[j, :])
            exprs = (λ[j].^U₀[j, :][nonzero_ids]).*vars[nonzero_ids]
            push!(free_action, collect(zip(vars[nonzero_ids], exprs)))
        end
    end
    mod_action = Vector{Tuple{Int, Vector{SparseAction}}}([])
    for (sᵢ, Uᵢ) in grading.mod_part
        sᵢ_action = Vector{SparseAction}([])
        if sᵢ == 2
            for j in axes(Uᵢ, 1)
                nonzero_ids = findall(!iszero, Uᵢ[j, :])
                push!(sᵢ_action, collect(zip(vars[nonzero_ids], -vars[nonzero_ids])))
            end
        else
            @var ω[sᵢ]
            for j in axes(Uᵢ, 1)
                nonzero_ids = findall(!iszero, Uᵢ[j, :])
                exprs = (ω[1].^Uᵢ[j, :][nonzero_ids]).*vars[nonzero_ids]
                push!(sᵢ_action, collect(zip(vars[nonzero_ids], exprs)))
            end
        end
        push!(mod_action, (sᵢ, sᵢ_action))
    end
    action = (free_action, mod_action)
    return ScalingGroup(grading, _structure(grading), vars, action)
end

Base.copy(s::ScalingGroup) = ScalingGroup(s.grading, s.structure, s.vars, s.action)

function Base.show(io::IO, scalings::ScalingGroup)
    action = scalings.action
    n_free, n_mod = length(action[1]), length(action[2])
    if n_free + n_mod == 0
        print(io, "ScalingGroup with 0 scalings")
        return
    end
    println(io, "ScalingGroup isomorphic to $(scalings.structure)")
    if n_free != 0
        print(io, " $(phrase(n_free, "free scaling")):")
        for free_action in scalings.action[1]
            print(io, "\n  ")
            for (j, (var, expr)) in enumerate(free_action)
                print(io, var, " ↦ ", expr)
                j < length(free_action) && print(io, ", ")
            end
        end
    end
    if n_mod != 0
        n_free != 0 && println(io, "\n")
        println(io, " modular scalings:")
        for (sᵢ, sᵢ_actions) in scalings.action[2]
            print(io, "  $(length(sᵢ_actions)) of order $(sᵢ):")
            for mod_action in sᵢ_actions
                print(io, "\n   ")
                for (j, (var, expr)) in enumerate(mod_action)
                    print(io, var, " ↦ ", expr)
                    j < length(mod_action) && print(io, ", ")
                end
            end
        end
    end
end

function _snf_scaling_symmetries(F::System)::Tuple{Vector{Int}, Vector{Matrix{Int}}}
    vars = vcat(F.variables, F.parameters)
    Es = [exponents_coefficients(f, vars)[1] for f in F.expressions]
    K = hcat([column_diffs(E) for E in Es]...)
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
    return s, Us
end

function _hnf_reduce(grading::Grading)
    U₀ = grading.free_part
    hnf_U₀ = isnothing(U₀) ? nothing : Matrix(hnf(matrix(ZZ, U₀)))
    hnf_Uᵢs = [(sᵢ, lift.(Matrix(hnf(matrix(GF(sᵢ), Uᵢ))))) for (sᵢ, Uᵢ) in grading.mod_part]
    return Grading(hnf_U₀, hnf_Uᵢs)
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
    vars = vcat(F.variables, F.parameters)
    s, U = _snf_scaling_symmetries(F)
    length(s) == 0 && return ScalingGroup(vars)
    grading = Grading(s, U)
    in_hnf && return ScalingGroup(_hnf_reduce(grading), vars)
    return ScalingGroup(grading, vars)
end

scaling_symmetries(
    F::SampledSystem;
    in_hnf::Bool=true
) = scaling_symmetries(F.system; in_hnf=in_hnf)

# TODO: extend to remove rows dependent on other blocks
function reduce(grading::Grading)
    hnf_grading = _hnf_reduce(grading)
    red_grading = Grading()
    U₀ = filter_rows(!iszero, hnf_grading.free_part)
    red_grading.free_part = size(U₀, 1) == 0 ? nothing : U₀
    for (sᵢ, Uᵢ) in hnf_grading.mod_part
        Uᵢ = filter_rows(!iszero, Uᵢ)
        size(Uᵢ, 1) != 0 && push!(red_grading.mod_part, (sᵢ, Uᵢ))
    end
end

function restrict_scalings(scalings::ScalingGroup, var_ids::Vector{Int})
    restr_grading = copy(scalings.grading)
    U₀ = restr_grading.free_part
    restr_grading.free_part = isnothing(U₀) ? nothing : U₀[:, var_ids]
    for (sᵢ, Uᵢ) in restr_grading.mod_part
        restr_grading.mod_part[i] = (sᵢ, Uᵢ[:, var_ids])
    end
    return ScalingGroup(reduce(restr_grading), scalings.vars[var_ids])
end

# TODO: what if vars is not a subset of scalings.vars?
function restrict_scalings(scalings::ScalingGroup, vars::Vector{Variable})
    var_ids = [findfirst(v->v==var, scalings.vars) for var in vars]
    if nothing in var_ids
        throw(ArgumentError("vars contains variables not present in scalings.vars"))
    end
    return restrict_scalings(scalings, var_ids)
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
    U₀ = grading.free_part
    deg_free = isnothing(U₀) ? Vector{Int}([]) : U₀*md
    return vcat(deg_free, [mod(Uᵢ*md, sᵢ) for (sᵢ, Uᵢ) in grading.mod_part]...)
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