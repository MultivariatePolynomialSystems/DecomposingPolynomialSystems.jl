using DecomposingPolynomialSystems
using Documenter

DocMeta.setdocmeta!(DecomposingPolynomialSystems, :DocTestSetup, :(using DecomposingPolynomialSystems); recursive=true)

makedocs(;
    modules=[DecomposingPolynomialSystems],
    authors="Viktor Korotynskiy <korotynskiy.viktor@gmail.com> and contributors",
    repo="https://github.com/azoviktor/DecomposingPolynomialSystems.jl/blob/{commit}{path}#{line}",
    sitename="DecomposingPolynomialSystems.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://azoviktor.github.io/DecomposingPolynomialSystems.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/azoviktor/DecomposingPolynomialSystems.jl",
    devbranch="main",
)
