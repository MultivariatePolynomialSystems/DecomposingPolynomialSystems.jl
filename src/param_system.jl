export ParametricSystem,
    unknowns,
    parameters,
    variables,
    nunknowns,
    nparameters,
    nvariables

struct ParametricSystem
    F::System
    unknowns::Vector{Variable}
    parameters::Vector{Variable}
end

function ParametricSystem(eqs::Vector{Expression};
                 unknowns::Vector{Variable}=Variable[],
                 parameters::Vector{Variable}=Variable[])
    if isempty(parameters)
        error("No parameters specified for ParametricSystem")
    end
    return ParametricSystem(
            System(eqs; variables=unknowns, parameters=parameters),
            unknowns,
            parameters
        )
end

unknowns(F::ParametricSystem) = F.unknowns
parameters(F::ParametricSystem) = F.parameters
variables(F::ParametricSystem) = vcat(F.unknowns, F.parameters)
nunknowns(F::ParametricSystem) = length(F.unknowns)
nparameters(F::ParametricSystem) = length(F.parameters)
nvariables(F::ParametricSystem) = nunknowns(F) + nparameters(F)

System(F::ParametricSystem) = F.F
equations(F::ParametricSystem) = expressions(System(F))
nequations(F::ParametricSystem) = length(System(F))

function Base.show(io::IO, F::ParametricSystem)
    println(io, "ParametricSystem with $(phrase(nequations(F), "polynomial"))")
    println(io, " $(phrase(nunknowns(F), "unknown")): ", join(unknowns(F), ", "))
    println(io, " $(phrase(nparameters(F), "parameter")): ", join(parameters(F), ", "))
    println(io, " $(phrase(nequations(F), "polynomial")):")
    for (i, eq) in enumerate(equations(F))
        print(io, "  $eq")
        i < nequations(F) && println(io)
    end
end

function run_monodromy(
    F::ParametricSystem,
    start_points::Union{Nothing, Tuple{AbstractVector{<:AbstractVector{<:Number}}, AbstractVector{<:Number}}}=nothing;
    options...
)
    if isnothing(start_points)
        MR = HC.monodromy_solve(System(F); permutations=true, options...)
    else
        sols, p₀ = start_points
        MR = HC.monodromy_solve(System(F), sols, ComplexF64.(p₀); permutations=true, options...)
    end
    if length(HC.solutions(MR)) == 1
        error("Only 1 solution found, no monodromy group available. Try running again...")
    end
    return SampledParametricSystem(F, MR)
end