export ExpressionMap

MiExpression = Union{Missing, Expression}

struct ExpressionMap
    domain_vars::Vector{Variable}
    image_vars::Vector{Variable}
    exprs::Vector{MiExpression}

    function ExpressionMap(domain_vars, image_vars, funcs)
        # TODO: exclude empty vectors, repetitions in vars
        return new(domain_vars, image_vars, funcs)
    end
end

# TODO: what if vars not in funcs? What if funcs has variables not present in vars?
function ExpressionMap(vars::Vector{Variable}, exprs::Vector{MiExpression})
    @assert length(vars) == length(exprs) "#vars ≂̸ #exprs, specify image variables"
    return ExpressionMap(vars, vars, exprs)
end

Base.getindex(f::ExpressionMap, i::Int) = (f.image_vars[i], f.exprs[i])
function Base.getindex(f::ExpressionMap, var::Variable)
    id = findfirst(x->x==var, f.image_vars)
    if isnothing(id)
        error("The variable $(var) isn't present in the image variables")
    end
    return f.exprs[id]
end

# TODO
function (f::ExpressionMap)(x)

end

# TODO
function Base.:(∘)(f::ExpressionMap, g::ExpressionMap)

end

# TODO
function is_dense(f::ExpressionMap)

end

function Base.show(io::IO, map::ExpressionMap)
    println(io, "ExpressionMap: ℂ$(superscriptnumber(length(map.domain_vars))) ⊃ X - - > ℂ$(superscriptnumber(length(map.exprs)))")
    println(io, " action:")
    if map.domain_vars == map.image_vars
        for (i, var) in enumerate(map.domain_vars)
            print(io, "  ", var, " ↦ ", map.exprs[i])
            i < length(map.domain_vars) && print(io, "\n")
        end
    else
        # TODO
    end
end