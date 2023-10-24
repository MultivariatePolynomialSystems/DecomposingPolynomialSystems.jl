using DecomposingPolynomialSystems

@var x a b
F = System([x^2+a*x+b]; variables=[x], parameters=[a,b])
F = run_monodromy(F)
mons = MonomialVector{Int8}([x,a,b], 3)
deck = symmetries_fixing_parameters!(F; degree_bound=1, param_dep=true)


using Test

@testset "DecomposingPolynomialSystems.jl" begin
    # Write your tests here.
end
