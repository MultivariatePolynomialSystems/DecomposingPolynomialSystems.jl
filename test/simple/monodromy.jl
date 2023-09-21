using DecomposingPolynomialSystems, HomotopyContinuation

@var x, y, p
F = System([x^2 + x + p, x + y + p]; parameters = [p])
xp0 = (x0, p0) = (ComplexF64.([1, 1]), ComplexF64.([-2]))
F = run_monodromy(F, xp0, max_loops_no_progress=3, tracker_options=TrackerOptions(parameters=CONSERVATIVE_TRACKER_PARAMETERS))
