using DecomposingPolynomialSystems

@var x, y, p
F = System([x^2 + x + p, x + y + p]; parameters = [p])
deck_transformations = compute_deck_transformations(F, degree=1);
deck_transformations[1]

@var x, a, b
F = System([x^2 + a*x + b]; parameters = [a,b])
# xp0 = (x0, p0) = (ComplexF64.([1]), ComplexF64.([1, -2]))
deck_transformations = compute_deck_transformations(F, degree=1);
deck_transformations[1]

@var x, a, b
F = System([x^2 + a*x + b]; parameters = [a, b])
xp0 = (ComplexF64.([1]), ComplexF64.([1, -2]))
F = run_monodromy(F, xp0)
deck_transformations = compute_deck_transformations!(F, degree=1, param_dep=true)

# PROBLEM. Try degree=1 and degree=2. On this variety 1/x and -x-p are equivalent.
@var x, p
F = System([x^2 + p*x + 1]; parameters = [p])
xp0 = (x0, p0) = (ComplexF64.([2]), ComplexF64.([-5/2]))
deck_transformations = compute_deck_transformations(F, xp0, degree=1, nd_perm=true);
deck_transformations[1]

@var x, y, p
F = System([x^2 + p*x + 1, x + y + p]; parameters = [p])
xp0 = (x0, p0) = (ComplexF64.([2, 1/2]), ComplexF64.([-5/2]))
deck_transformations = compute_deck_transformations(F, xp0, degree=1, nd_perm=false);
deck_transformations[1]


# PROBLEM. Try degree=1 and degree=2. On this variety, a/x and -1-x are equivalent.
@var x, a
F = System([x^2 + x + a]; parameters = [a])
xp0 = (x0, p0) = (ComplexF64.([1]), ComplexF64.([-2]))
deck_transformations = compute_deck_transformations(F, xp0, degree=1);
deck_transformations[1]

@var x, a, b
F = System([x^4 + a*x^2 + b]; parameters = [a,b])
xp0 = (x0, p0) = (ComplexF64.([1]), ComplexF64.([1, -2]))
deck_transformations = compute_deck_transformations(F, xp0, degree=3);
deck_transformations[1]

@var x y a b
F = System([x^4 + y - 2*a,  x^4 + x^2 - 2*b*y^2]; parameters = [a,b])
xp0 = (ComplexF64.([1, 1]), ComplexF64.([1, 1]))
deck_transformations = compute_deck_transformations(F, xp0, degree=3);
deck_transformations[1]

@var x y p[1:2]
F = System([x^4 + y^2 - 2p[1],  x^4 + x^2 - 2*p[2]*y^2]; parameters = p)
xp0 = (ComplexF64.([1, 1]), ComplexF64.([1, 1]))
deck_transformations = compute_deck_transformations(F, xp0, degree=3);
deck_transformations[1]

@var x y p[1:2]
F = System([x^2 - y^2 - p[1], 2*x*y - p[2]]; parameters = p)
xp0 = (ComplexF64.([1, 0]), ComplexF64.([1, 0]))
deck_transformations = compute_deck_transformations(F, xp0, degree=3);
deck_transformations[1]

# Numerical Problems for degree > 5
@var x p
F = System([2*x^4 + 3*x^3 + 1/p*x^2 + 3*x + 2]; parameters = [p])
xp0 = (ComplexF64.([-2]), ComplexF64.([-1]))
deck_transformations = compute_deck_transformations(F, xp0, degree=5);
deck_transformations[1]

@var x p
F = System([2*(x+1)^4 + 3*(x+1)^3 + 1/p*(x+1)^2 + 3*(x+1) + 2]; parameters = [p])
xp0 = (ComplexF64.([-3]), ComplexF64.([-1]))
deck_transformations = compute_deck_transformations(F, xp0, degree=3);
deck_transformations[1]

# Here x[2] is fixed to 1, causes problems
@var x[1:2] p
F = System([x[1]^4*x[2] + x[1]^2*x[2]^3 + p*x[2]^5, x[2]-1]; parameters = [p])
xp0 = (ComplexF64.([1, 1]), ComplexF64.([-2]))
deck_transformations = compute_deck_transformations(F, xp0, degree=2);
deck_transformations[1]