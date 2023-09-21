include("../../../../src/implicitization/implicitize.jl")
include("../eqs_wo_depths.jl")

F = eqs_5pp()
xz0 = (x0, z0) = fabricateSolution()
norm(F(x0, z0))

F = run_monodromy(F, xz0)
F = sample_system(F, 100)
