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
        collapselevel=2
    ),
    pages = [
        "Introduction" => "index.md",
        "Sampling Polynomial Systems" => "sampling.md",
        "Symmetries" => [
            "Scaling symmetries" => "symmetries/scalings.md",
            "Deck transformations" => "symmetries/deck.md",
        ],
        "Invariants" => [
            "Invariants" => "invariants/invariants.md",
        ],
        "Decomposing Polynomial Systems" => [
            "Decompose" => "decomposition/decompose.md",
        ]
    ],
)

deploydocs(;
    repo="github.com/MultivariatePolynomialSystems/DecomposingPolynomialSystems.jl.git",
    devbranch="main",
)
