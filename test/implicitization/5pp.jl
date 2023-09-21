include("../../../test/computer_vision/5pp/eqs_wo_depths.jl")

F = eqs_5pp()
(x0, p0) = fabricateSolution()
norm(F(x0, p0))

unknowns = variables(F)
params = parameters(F)
