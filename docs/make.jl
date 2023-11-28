push!(LOAD_PATH,"../src/")

using DecomposingPolynomialSystems
using Documenter

DocMeta.setdocmeta!(DecomposingPolynomialSystems, :DocTestSetup, :(using DecomposingPolynomialSystems); recursive=true)

makedocs(;
    modules=[DecomposingPolynomialSystems],
    authors="Viktor Korotynskiy <korotynskiy.viktor@gmail.com> and contributors",
    repo="https://github.com/MultivariatePolynomialSystems/DecomposingPolynomialSystems.jl/blob/{commit}{path}#{line}",
    sitename="DecomposingPolynomialSystems.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://multivariatepolynomialsystems.github.io/DecomposingPolynomialSystems.jl",
        edit_link="main",
        assets=["assets/custom.css"],
        collapselevel=1
    ),
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
            "Monomials" => "DPS/monomials.md",
            "Symmetries" => "DPS/symmetries.md",
            "Invariants" => "DPS/invariants.md",
        ],
        "Examples from Algebraic Vision" => [
            "5-point problem" => "AV/5pp.md",
        ],
    ],
)

deploydocs(;
    repo="github.com/MultivariatePolynomialSystems/DecomposingPolynomialSystems.jl.git",
    devbranch="main",
)
