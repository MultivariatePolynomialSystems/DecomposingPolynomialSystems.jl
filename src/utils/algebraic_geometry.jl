using HomotopyContinuation: differentiate
export jac

function jac(F::System)
    return differentiate(F.expressions, F.variables)
end

function jac(F::System, x0::Vector{CC})
    return Matrix{CC}(subs(differentiate(F.expressions, F.variables), F.variables => x0))
end

function jac(F::System, xp0::Tuple{Vector{CC}, Vector{CC}})
    return Matrix{CC}(subs(differentiate(F.expressions, F.variables), vcat(F.variables, F.parameters) => vcat(xp0[1], xp0[2])))
end