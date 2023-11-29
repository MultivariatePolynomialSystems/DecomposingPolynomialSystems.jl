export Grading,
    ScalingGroup,
    scaling_symmetries

using AbstractAlgebra: ZZ, matrix, GF, lift, hnf, snf_with_transform 
using LinearAlgebra: diag

mutable struct Grading{Tv<:Integer,Ti<:Integer}
    nscalings::Int
    free_part::Union{Nothing, SparseMatrixCSC{Tv,Ti}}
    mod_part::Vector{Tuple{Tv, SparseMatrixCSC{Tv,Ti}}}
end

Grading{Tv,Ti}() where {Tv<:Integer,Ti<:Integer} = Grading{Tv,Ti}(0, nothing, [])

function Grading{Tv,Ti}(
    s::Vector{Int},
    U::Vector{Matrix{Int}}
) where {Tv<:Integer,Ti<:Integer}

    grading = Grading{Tv,Ti}()
    for (sᵢ, Uᵢ) in zip(s, U)
        if sᵢ == 0
            grading.free_part = Uᵢ
        else
            push!(grading.mod_part, (sᵢ, Uᵢ))
        end
        grading.nscalings += size(Uᵢ, 1)
    end
    return grading
end

nfree(grading::Grading) = isnothing(grading.free_part) ? 0 : size(grading.free_part, 1)
nscalings(grading::Grading) = grading.nscalings
Base.isempty(grading::Grading) = isnothing(grading.free_part) && isempty(grading.mod_part)
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

"""
    ScalingGroup{Tv<:Integer,Ti<:Integer}

A `ScalingGroup` is the result of `scaling_symmetries` computation.
"""
struct ScalingGroup{Tv<:Integer,Ti<:Integer}
    grading::Grading{Tv,Ti}
    structure::String
    vars::Vector{Variable}
    action::Tuple{Vector{SparseAction}, Vector{Tuple{Tv, Vector{SparseAction}}}}
end

ScalingGroup{Tv,Ti}(
    vars::Vector{Variable}
) where {Tv<:Integer,Ti<:Integer} = ScalingGroup(Grading{Tv,Ti}(), "N/A", vars, ([], []))

function ScalingGroup(
    grading::Grading{Tv,Ti},
    vars::Vector{Variable}
) where {Tv<:Integer,Ti<:Integer}
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
            nzind, nzval = findnz(U₀[j, :])
            exprs = (λ[j].^nzval).*vars[nzind]
            push!(free_action, collect(zip(vars[nzind], exprs)))
        end
    end
    mod_action = Vector{Tuple{Tv, Vector{SparseAction}}}([])
    for (sᵢ, Uᵢ) in grading.mod_part
        sᵢ_action = Vector{SparseAction}([])
        if sᵢ == 2
            for j in axes(Uᵢ, 1)
                nzind, _ = findnz(Uᵢ[j, :])
                push!(sᵢ_action, collect(zip(vars[nzind], -vars[nzind])))
            end
        else
            @var ω[Int(sᵢ)]
            for j in axes(Uᵢ, 1)
                nzind, nzval = findnz(Uᵢ[j, :])
                exprs = (ω[1].^nzval).*vars[nzind]
                push!(sᵢ_action, collect(zip(vars[nzind], exprs)))
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
        for (i, (sᵢ, sᵢ_actions)) in enumerate(scalings.action[2])
            print(io, "  $(length(sᵢ_actions)) of order $(sᵢ):")
            for mod_action in sᵢ_actions
                print(io, "\n   ")
                for (j, (var, expr)) in enumerate(mod_action)
                    print(io, var, " ↦ ", expr)
                    j < length(mod_action) && print(io, ", ")
                end
            end
            i < length(scalings.action[2]) && print(io, "\n")
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

function _hnf_reduce!(
    s::Vector{Int},
    U::Vector{Matrix{Int}}
)
    for (i, (sᵢ, Uᵢ)) in enumerate(zip(s, U))
        if sᵢ == 0
            U[i] = Matrix(hnf(matrix(ZZ, Uᵢ)))
        else
            try
                U[i] = Int.(lift.(Matrix(hnf(matrix(GF(sᵢ), Uᵢ)))))
            catch
                U[i] = mod.(Uᵢ, sᵢ)
            end
        end
    end
    return s, U
end

function _hnf_reduce(grading::Grading{Tv,Ti}) where {Tv<:Integer,Ti<:Integer}
    U₀ = grading.free_part
    red_grading = Grading{Tv,Ti}()
    if !isnothing(U₀)
        red_grading.free_part = Matrix(hnf(matrix(ZZ, U₀)))
        red_grading.nscalings = size(U₀, 1)
    end
    for (sᵢ, Uᵢ) in grading.mod_part
        try
            Uᵢ = lift.(Matrix(hnf(matrix(GF(sᵢ), Uᵢ))))
        catch
            Uᵢ = mod.(Uᵢ, sᵢ)
        end
        push!(red_grading.mod_part, (sᵢ, Uᵢ))
        red_grading.nscalings += size(Uᵢ, 1)
    end
    return red_grading
end

"""
    scaling_symmetries(F::System)

Given a polynomial system `F` returns the group of scaling symmetries 
of `F`. The scalings that change the parameters are considered as well.

```julia-repl
julia> @var x y a b c;

julia> F = System([x^4+a^2+1, y^2+b+c]; variables=[x, y], parameters=[a,b,c]);

julia> scaling_symmetries(F)
ScalingGroup isomorphic to ℤ × ℤ₄ × ℤ₂
 1 free scaling:
  y ↦ y*λ, b ↦ b*λ^2, c ↦ c*λ^2

 modular scalings:
  1 of order 4:
   x ↦ x*ω₄^3, y ↦ y*ω₄, b ↦ b*ω₄^2, c ↦ c*ω₄^2
  1 of order 2:
   x ↦ -x, y ↦ -y, a ↦ -a
```
"""
function scaling_symmetries(F::System)
    vars = vcat(F.variables, F.parameters)
    s, U = _snf_scaling_symmetries(F)
    length(s) == 0 && return ScalingGroup{Int8, Int16}(vars)
    _hnf_reduce!(s, U)
    return ScalingGroup(Grading{Int8, Int16}(s, U), vars)
end

scaling_symmetries(F::SampledSystem) = scaling_symmetries(F.system)

# TODO: extend to remove rows dependent on other blocks
function reduce(grading::Grading{Tv,Ti}) where {Tv<:Integer,Ti<:Integer}
    hnf_grading = _hnf_reduce(grading)
    red_grading = Grading{Tv,Ti}()
    U₀ = take_rows(!iszero, hnf_grading.free_part)
    if size(U₀, 1) != 0
        red_grading.free_part = U₀
        red_grading.nscalings = size(U₀, 1)
    end
    for (sᵢ, Uᵢ) in hnf_grading.mod_part
        Uᵢ = take_rows(!iszero, Uᵢ)
        if size(Uᵢ, 1) != 0
            push!(red_grading.mod_part, (sᵢ, Uᵢ))
            red_grading.nscalings += size(Uᵢ, 1)
        end
    end
    return red_grading
end

function restrict_scalings(scalings::ScalingGroup, var_ids::Vector{Int})
    restr_grading = copy(scalings.grading)
    U₀ = restr_grading.free_part
    restr_grading.free_part = isnothing(U₀) ? nothing : U₀[:, var_ids]
    for (i, (sᵢ, Uᵢ)) in enumerate(restr_grading.mod_part)
        restr_grading.mod_part[i] = (sᵢ, Uᵢ[:, var_ids])
    end
    return ScalingGroup(reduce(restr_grading), scalings.vars[var_ids])
end

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

function HC.degree(
    mexp::SparseVector{Tv,Ti},
    grading::Grading{Tv,Ti}
) where {Tv<:Integer,Ti<:Integer}
    deg = spzeros(Tv, Ti, nscalings(grading))
    U₀ = grading.free_part
    if !isnothing(U₀)
        deg[1:size(U₀,1)] = U₀*mexp
    end
    k = nfree(grading)
    for (sᵢ, Uᵢ) in grading.mod_part
        deg[(k+1):(k+size(Uᵢ,1))] = mod.(Uᵢ*mexp, sᵢ)
        k += size(Uᵢ, 1)
    end
    return deg
end

function to_classes(
    mons::DenseMonomialVector{Tv,Ti},
    grading::Grading{Tv,Ti}
) where {Tv<:Integer,Ti<:Integer}

    classes = Dict{SparseVector{Tv,Ti}, SparseMonomialVector{Tv,Ti}}()
    for (uexp, pexp) in mons
        deg = HC.degree(vcat(uexp, pexp), grading)
        if isnothing(get(classes, deg, nothing)) # the key doesn't exist
            classes[deg] = SparseMonomialVector{Tv,Ti}(mons.unknowns, mons.parameters)
        end
        push!(classes[deg], (uexp, pexp))
    end
    return classes
end
