using DecomposingPolynomialSystems
using Test


# @var x, a, b
# F = System([x^3 + a*x + b]; parameters = [a,b])
# xz0 = (x0, z0) = (ComplexF64.([2]), ComplexF64.([1, -10]))
# F = decompose(F, xz0)
# @test F.is_decomposable == false


# @var x, a, b
# F = System([a*x^4 + b*x^2 + a]; parameters = [a,b])
# xz0 = (x0, z0) = (ComplexF64.([2]), ComplexF64.([1, -17/4]))
# (F, factorizing_maps) = compute_factorizing_maps(F, xz0, degree=2)
# F = sample_system(F, 10)
# implicitize(F, factorizing_maps[3], F.block_partitions[3], degree=3)
# F = decompose(F, xz0, degDeck=1, degFact=2, degImpl=2)
# @test F.is_decomposable == true


@var x, a, b
F = System([x^4 + a*x^2 + b]; parameters = [a,b])
xz0 = (x0, z0) = (ComplexF64.([2]), ComplexF64.([1, -20]))
F = decompose(F, xz0, degDeck=1, degFact=2, degImpl=2)
@test F.is_decomposable == true
