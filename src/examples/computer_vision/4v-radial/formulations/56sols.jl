include("../../../../utils/utils.jl")
using DecomposingPolynomialSystems
using HomotopyContinuation, BlockDiagonals

function radial_uncalibrated_eqs()
    @var p[1:4, 1:4], l[1:2,1:13,1:4]
    id4 = eye(4)
    Ps = vcat([[id4[1,:]'; [p[1,1] p[1,1] p[1,1] p[1,1]]]], [[id4[i,:]'; [p[i,j] for j=1:4]'] for i=2:4])
    stackedPs = Matrix(hcat(Ps'...)')
    Ms = [hcat(stackedPs, BlockDiagonal(M2VM(l[:,i,:]))) for i=1:13]
    eqs = [exprDet(Ms[i]) for i = 1:13]
    return System(eqs; variables=vcat([p[1,1]], reshape(p[2:4,:], 12)), parameters=reshape(l, 2*13*4))
end

function fabricateProblemSolutionPair()
    id4 = eye(4)
    p = transpose(hcat(randn(ComplexF64)*ones(4), randn(ComplexF64, 4), randn(ComplexF64, 4), randn(ComplexF64, 4)))
    Ps = [[id4[i,:]'; transpose(p[i,:])] for i=1:4]
    stackedPs = Matrix(transpose(hcat(transpose(Ps)...)))
    X = rand(ComplexF64, 3, 13)
    l = hcat([Ps[i]*a2p(X) for i=1:4]...)
    return (vcat([p[1,1]], reshape(p[2:4,:], 12)), reshape(l, 2*13*4))
end

F = radial_uncalibrated_eqs()
(x0, p0) = fabricateProblemSolutionPair()
norm(F(x0,p0))


F = run_monodromy(F, (x0,p0))
vars = variables(F.equations)
p11,p21,p31,p41,p22,p32,p42,p23,p33,p43,p24,p34,p44 = vars
A = [1,p11,p22,1/p33,1/p44,p21,p31,1/p41,
     p23*p32,p24*p42,1/(p34*p43),p23*p31+p21*p32,1/(p24*p41+p21*p42),p34*p41+p31*p43,

     p21^2, p31^2, p41^2, p11^3, p21^3, p31^3, p41^3, p21*p31*p41,
     
     1/p22, 1/p21, 1/p31, 1/p21^2, 1/p31^2, 1/p41^2]
are_LI(F, A)
B = [p23, p23*p24, p23*p33, p23*p34, p23*p41, p23*p42, p23*p43, p23*p44,
     p24,p24^2,p24*p31,p24*p32,p24*p33,p24*p34,p24*p43,p24*p44,p32,p32^2,
     p32*p33,p32*p34,p32*p41,p32*p42,p32*p43,p32*p44,p34,p34^2,p34*p42,p34*p44]
Binv = vcat(A, B)
are_LI(F, Binv)

eval_at_sols(F, [p23*p31+p21*p32])
eval_at_sols(F, [p34*p41+p31*p43])
eval_at_sols(F, [p24*p41+p21*p42])


monodromy_permutations = filter_permutations_with_zeros(HomotopyContinuation.permutations(MR))
G = permutations_to_group(monodromy_permutations)
GAP.Globals.Order(G)

A56 = GAP.Globals.AlternatingGroup(56)
GAP.Globals.IsSubgroup(A56,G)
