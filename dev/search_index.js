var documenterSearchIndex = {"docs":
[{"location":"SPS/symmetries/#Symmetries","page":"Symmetries","title":"Symmetries","text":"","category":"section"},{"location":"SPS/symmetries/","page":"Symmetries","title":"Symmetries","text":"Some decompositions are caused by the presence of non-trivial symmetries of the polynomial system that are given by rational maps. There are 2 important types of symmetries. The 1st type is scaling symmetries, namely, the maps that act by scaling individual variables of the formuation. The 2nd type are deck transformations, namely, the birational maps that fix the parameters of the polynomial system.","category":"page"},{"location":"SPS/symmetries/#Scaling-symmetries","page":"Symmetries","title":"Scaling symmetries","text":"","category":"section"},{"location":"SPS/symmetries/","page":"Symmetries","title":"Symmetries","text":"scaling_symmetries","category":"page"},{"location":"SPS/symmetries/#DecomposingPolynomialSystems.scaling_symmetries","page":"Symmetries","title":"DecomposingPolynomialSystems.scaling_symmetries","text":"scaling_symmetries(F::System; in_hnf::Bool=true)\n\nGiven a polynomial system F returns the group of scaling symmetries  of the polynomial system F.\n\njulia> @var x[1:2] p[1:2];\n\njulia> F = System([x[1]^2 - x[2]^2 - p[1], 2*x[1]*x[2] - p[2]]; variables=x, parameters=p);\n\njulia> scalings = scaling_symmetries(F)\nScalingSymmetryGroup with 2 scalings\n infinite scalings: 1\n finite scalings:\n  1 of order 2\n vars: x₁, x₂, p₁, p₂\n\n\n\n\n\n","category":"function"},{"location":"SPS/symmetries/#Deck-transformations","page":"Symmetries","title":"Deck transformations","text":"","category":"section"},{"location":"SPS/symmetries/","page":"Symmetries","title":"Symmetries","text":"symmetries_fixing_parameters","category":"page"},{"location":"SPS/symmetries/#DecomposingPolynomialSystems.symmetries_fixing_parameters","page":"Symmetries","title":"DecomposingPolynomialSystems.symmetries_fixing_parameters","text":"symmetries_fixing_parameters(F::System; degree_bound=1, param_dep=true, tol=1e-5)\n\nGiven a polynomial system F returns the group of symmetries  of the polynomial system F that fix the parameters. The keyword argument degree_bound is used to set the upper bound for the degrees of numerator and denominator polynomials in expressions for the symmetries.\n\njulia> @var x[1:2] p[1:2];\n\njulia> F = System([x[1]^2 - x[2]^2 - p[1], 2*x[1]*x[2] - p[2]]; variables=x, parameters=p);\n\njulia> symmetries_fixing_parameters(F; degree_bound=1, param_dep=false)\nDeckTransformationGroup of order 4\n structure: C2 x C2\n action:\n  1st map:\n   x₁ ↦ x₁\n   x₂ ↦ x₂\n  2nd map:\n   x₁ ↦ (0.0 + 1.0*im)*x₂\n   x₂ ↦ (0.0 - 1.0*im)*x₁\n  3rd map:\n   x₁ ↦ (-1.0 + 0.0*im)*x₁\n   x₂ ↦ (-1.0 + 0.0*im)*x₂\n  4th map:\n   x₁ ↦ (0.0 - 1.0*im)*x₂\n   x₂ ↦ (0.0 + 1.0*im)*x₁\n\n\n\n\n\n","category":"function"},{"location":"AV/5pp/#Problem-Description","page":"5-point problem","title":"Problem Description","text":"","category":"section"},{"location":"#Introduction","page":"Introduction","title":"Introduction","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"DecomposingPolynomialSystems.jl is a Julia package for decomposing systems of polynomial equations, i.e. representing a possibly complicated polynomial system as a sequence of simpler polynomial subsystems.","category":"page"},{"location":"#Package-Features","page":"Introduction","title":"Package Features","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Computing symmetries that fix the parameters of the polynomial system\nComputing invariants... ","category":"page"},{"location":"NAG/symmetries/#Symmetries","page":"Symmetries","title":"Symmetries","text":"","category":"section"},{"location":"NAG/monodromy/#Monodromy","page":"Monodromy","title":"Monodromy","text":"","category":"section"}]
}
