include("../../../src/utils/utils.jl")
using HomotopyContinuation, BlockDiagonals

function metric_upgrade_eqs()
    @var p[1:4, 1:4], q[1:10], l[1:2,1:13,1:4]
    id4 = eye(4)
    Ps = vcat([[id4[1,:]'; [p[1,1] p[1,1] p[1,1] p[1,1]]]], [[id4[i,:]'; [p[i,j] for j=1:4]'] for i=2:4])
    stackedPs = Matrix(hcat(Ps'...)')
    Ms = [hcat(stackedPs, BlockDiagonal(M2VM(l[:,i,:]))) for i=1:13]
    eqs1 = [exprDet(Ms[i]) for i = 1:13]
    Q = [[q[1] q[2] q[3] q[4]]; [q[2] q[5] q[6] q[7]]; [q[3] q[6] q[8] q[9]]; [q[4] q[7] q[9] q[10]]]
    eqs2 = [M2V(Ps[i]*Q*transpose(Ps[i])) for i=1:4]
    new_eqs2 = []
    for eqs in eqs2
        append!(new_eqs2, [eqs[1]-eqs[4], eqs[2]])
    end
    eqs = vcat(vcat(eqs1, new_eqs2), [exprDet(Q), q[1]-1])
    return System(eqs; variables=vcat([p[1,1]], reshape(p[2:4,:], 12), q), parameters=reshape(l, 2*13*4))
end

function fabricateProblemSolutionPair()
    cs = randn(ComplexF64, 3, 4)
    ts = randn(ComplexF64, 2, 4)
    Ps = [ct2sPrad(cs[:,i], ts[:,i]) for i=1:4]

    A = vcat([transpose(Ps[i][1,:]) for i=1:4]...)
    d = zeros(ComplexF64, 4)
    for i = 1:4
        As = copy(A)
        As[i,:] = 1/Ps[1][1,1]*Ps[1][2,:]
        d[i] = det(As)
    end
    d = d/norm(d)
    H = inv(diagm(d)*A)
    Ps = [Ps[i]*H for i=1:4]
    Q = inv(H)*diagm([1,1,1,0])*inv(transpose(H))

    Ps = [1/Ps[i][1,i]*Ps[i] for i=1:4]

    s = zeros(ComplexF64, 4)
    for i = 1:4
        sI = Ps[i]*Q*transpose(Ps[i])
        s[i] = sI[1,1]
    end

    X = rand(ComplexF64, 3, 13)
    l = hcat([Ps[i]*a2p(X) for i=1:4]...)
    P = vcat([transpose(Ps[i][2,:]) for i=1:4]...)

    q = [Q[1,1], Q[1,2], Q[1,3], Q[1,4], Q[2,2], Q[2,3], Q[2,4], Q[3,3], Q[3,4], Q[4,4]]
    q = 1/q[1]*q

    return (vcat([P[1,1]], reshape(P[2:4,:], 12), q), reshape(l, 2*13*4))
end

F = metric_upgrade_eqs()
F.expressions
(x0, p0) = fabricateProblemSolutionPair()
norm(F(x0,p0))

MR = monodromy_solve(F, [x0], p0, permutations=true, max_loops_no_progress=10; tracker_options=TrackerOptions(parameters=:fast))

# corrector_steps = 1
# min_step_size
# max_loops_no_progress --> higher

ps = permutations(MR)
monodromy_permutations = filter_permutations_with_zeros(HomotopyContinuation.permutations(MR))
G = permutations_to_group(monodromy_permutations)
GAP.Globals.Order(G)
