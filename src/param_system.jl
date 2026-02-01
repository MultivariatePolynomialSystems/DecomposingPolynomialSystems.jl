export ParametricSystem,
    unknowns,
    parameters,
    variables,
    nunknowns,
    nparameters,
    nvariables

const MONODROMY_SOLVE_REF = "https://www.juliahomotopycontinuation.org/HomotopyContinuation.jl/stable/monodromy/"

struct ParametricSystem
    system::System
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

System(F::ParametricSystem) = F.system
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

(F::ParametricSystem)(
    x₀::AbstractVector{<:Number},
    p₀::AbstractVector{<:Number}
) = F.system(x₀, p₀)


"""
    run_monodromy(F::ParametricSystem, start_points=nothing; options...) -> SampledParametricSystem

Runs [`monodromy_solve`]($(MONODROMY_SOLVE_REF)) on a given polynomial system `F` with starting
solutions `start_points[1]` and parameters `start_points[2]` (if given).

```julia-repl
julia> @var x a b;

julia> F = ParametricSystem([x^3+a*x+b]; unknowns=[x], parameters=[a,b]);

julia> F = run_monodromy(F, ([[1]], [1,-2]); max_loops_no_progress = 10)
SampledParametricSystem with 3 samples
 1 unknown: x
 2 parameters: a, b

 number of solutions: 3
 sampled instances: 1
```
"""
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

function track_parameter_homotopy(
    F::ParametricSystem,
    (x₀, p₀)::NTuple{2, AbstractVector{<:Number}},
    p₁::AbstractVector{<:Number},
    p_inter::AbstractVector{<:Number} # intermediate parameter
)

    H₁ = ParameterHomotopy(System(F); start_parameters=p₀, target_parameters=p_inter)
    res = track(Tracker(H₁), x₀)
    if !is_success(res)
        @warn "Tracking was not successful: stopped at t = $(res.t)"
    end

    H₂ = ParameterHomotopy(System(F); start_parameters=p_inter, target_parameters=p₁)
    res = track(Tracker(H₂), solution(res))
    if !is_success(res)
        @warn "Tracking was not successful: stopped at t = $(res.t)"
    end

    return solution(res)
end