using DecomposingPolynomialSystems

@var x, a, b
F = System([x^3 + a*x + b]; parameters = [a,b])
x0 = ComplexF64.([2])
p0 = ComplexF64.([1, -10])
xp0 = (x0, p0)


@var x, a, b
F = System([a*x^4 + b*x^2 + a]; parameters = [a,b])
x0 = ComplexF64.([2])
p0 = ComplexF64.([-4/17, 1])
xp0 = (x0, p0)
compute_factorizing_maps(F, xp0, degree=2)


@var x, a, b
F = System([a*x^8 + b*x^4 + a]; parameters = [a,b])
x0 = ComplexF64.([2])
p0 = ComplexF64.([-16/257, 1])
xp0 = (x0, p0)
compute_deck_transformations(F, xp0, degree=2)
compute_factorizing_maps(F, xp0, degree=4)