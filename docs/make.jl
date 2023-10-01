push!(LOAD_PATH,"../src/")
using DecomposingPolynomialSystems, Documenter
makedocs(
    sitename = "DecomposingPolynomialSystems.jl",
    modules  = [DecomposingPolynomialSystems],
    format=Documenter.HTML(; collapselevel=1),
    pages = [
        "Introduction" => "index.md",
        "Numerical Algebraic Geometry" => [
            "Monodromy" => "NAG/monodromy.md",
            "Symmetries" => "NAG/symmetries.md",
        ],
        "Data Types" => [
            "SampledSystem" => "Types/sampled_system.md",
            "SymmetryGroup" => "Types/symmetry_group.md",
        ],
        "Decomposing Polynomial Systems" => [
            "Symmetries" => "SPS/symmetries.md",
            "Invariants" => "SPS/invariants.md",
        ],
        "Examples from Algebraic Vision" => [
            "5-point problem" => "AV/5pp.md",
        ],
    ])
deploydocs(;
    repo="github.com/MultivariatePolynomialSystems/DecomposingPolynomialSystems.jl.git",
)