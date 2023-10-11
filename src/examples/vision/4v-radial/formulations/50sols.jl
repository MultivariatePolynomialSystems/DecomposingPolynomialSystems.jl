include("../../../src/utils/utils.jl")
using HomotopyContinuation, BlockDiagonals

function radial_7pt_eqs()
    @var y[1:3], t[1:2,3:4], l[1:2,1:7,1:4]
    cs = vcat([[0,0,0]], [[0, y[i], 0] for i=1:3])
    ts = [[0;0], [0;1], t[:,1], t[:,2]]
    Ps = [ct2sPrad(cs[i], ts[i]) for i = 1:4]
    stackedPs = Matrix(hcat(Ps'...)')
    Ms = [hcat(stackedPs, BlockDiagonal(M2VM(l[:,i,:]))) for i=1:7]
    eqs = [exprDet(Ms[i]) for i = 1:7]
    return System(eqs; variables=vcat(y, reshape(t,4)), parameters=reshape(l, 2*7*4))
end

function fabricateProblemSolutionPair()
    y = rand(ComplexF64, 3)
    cs = vcat([[0,0,0]], [[0, y[i], 0] for i=1:3])
    ts = hcat([0;0], [0;1], rand(ComplexF64, 2, 2))
    Ps = [ct2sPrad(cs[i], ts[:,i]) for i=1:4]
    X = rand(ComplexF64, 3, 7)
    l = hcat([Ps[i]*a2p(X) for i=1:4]...)
    return (vcat(y, reshape(ts[:,3:4], 4)), reshape(l, 2*7*4))
end

(x0, p0) = fabricateProblemSolutionPair()
F = radial_7pt_eqs()
norm(F(x0,p0))


J = subs(differentiate(F.expressions, F.variables), vcat(F.variables, F.parameters) => vcat(x0, p0))
J = Matrix{ComplexF64}(J)
det(J)

MR = monodromy_solve(F, [x0], p0, permutations=true, max_loops_no_progress=5)
solutions(MR)

monodromy_permutations = filter_permutations_with_zeros(HomotopyContinuation.permutations(MR))
G = permutations_to_group(monodromy_permutations)
GAP.Globals.StructureDescription(G)
GAP.Globals.Order(G)
