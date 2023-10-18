include("../../../../utils/utils.jl")
using HomotopyContinuation, BlockDiagonals

function radial_13pt_eqs()
    @var x[2:4], y[2:4], z[2:4], t[1:2,3:4], l[1:2,1:13,1:4]
    cs = vcat([[0,0,0]], [[x[i], y[i], z[i]] for i=1:3])
    ts = [[0;0], [1;0], t[:,1], t[:,2]]
    Ps = [ct2sPrad(cs[i], ts[i]) for i = 1:4]
    stackedPs = Matrix(hcat(Ps'...)')
    Ms = [hcat(stackedPs, BlockDiagonal(M2VM(l[:,i,:]))) for i=1:13]
    eqs = [exprDet(Ms[i]) for i = 1:13]
    return System(eqs; variables=vcat(x, y, z, reshape(t,4)), parameters=reshape(l, 2*13*4))
end

function fabricateProblemSolutionPair()
    x, y, z = M2VV(randn(CC, 3, 3))
    ts = [[0;0], [1;0], randn(CC, 2), randn(CC, 2)]
    Rs = vcat([Matrix(I(3))], [c2R([x[i]; y[i]; z[i]]) for i=1:3])
    Ps = [hcat(Rs[i][1:2,:], ts[:,i]) for i = 1:4]
    Ps = [ct2sPrad(cs[i], ts[i]) for i = 1:4]
    stackedPs = Matrix(hcat(Ps'...)')
    X = rand(ComplexF64, 3, 13)
    l = hcat([Ps[i]*a2p(X) for i=1:4]...)
    Ms = [hcat(stackedPs, BlockDiagonal([V2M(l[:,i]), V2M(l[:,13+i]), V2M(l[:,2*13+i]), V2M(l[:,3*13+i])])) for i=1:13]
    return (vcat(x, y, z, reshape(ts[:,3:4], 4)), reshape(l, 2*13*4))
end

(x0, p0) = fabricateProblemSolutionPair()
F = radial_13pt_eqs()
norm(F(x0,p0))


J = subs(differentiate(F.expressions, F.variables), vcat(F.variables, F.parameters) => vcat(x0, p0))
J = Matrix{ComplexF64}(J)
det(J)

MR = monodromy_solve(F, [x0], p0, permutations=true, max_loops_no_progress=5)
# sols = solutions(MR)
# params = parameters(MR)

using FileIO, JLD2
sols = FileIO.load("test/computer_vision/calibrated/4v13p_radial/4v13p_radial_solutions.jld2", "sols")
params = FileIO.load("test/computer_vision/calibrated/4v13p_radial/4v13p_radial_parameters.jld2", "params")

function evaluateCameras(sol)
    ts = hcat([0;0], [1;0], sol[10:11], sol[12:13])
    Rs = vcat([Matrix(I(3))], [cayley_rot(sol[((0:2)*3).+i]) for i=1:3])
    return [hcat(Rs[i][1:2,:], ts[:,i]) for i = 1:4]
end

function evaluateRQT(sol)
    ts = hcat([0;0], [1;0], sol[10:11], sol[12:13])
    Rs = vcat([Matrix(I(3))], [cayley_rot(sol[((0:2)*3).+i]) for i=1:3])
    Ps = [hcat(Rs[i][1:2,:], ts[:,i]) for i = 1:4]
    return [(-1)^(i1+i2+i3+i4)*det([transpose(Ps[1][i1,:]); transpose(Ps[2][i2,:]); transpose(Ps[3][i3,:]); transpose(Ps[4][i4,:])]) for (i1,i2,i3,i4) in vec(reshape(collect(Iterators.product(ntuple(_ -> 1:2, 4)...)), 2^4, 1))]
end

function larssonSymmetry(sol, camera_id)
    new_sol = copy(sol)
    xi,yi,zi = sol[((0:2)*3).+(camera_id-1)]
    new_sol[((0:2)*3).+(camera_id-1)] = [-yi/zi, xi/zi, -1/zi] # my [x[i]; y[i]; z[i]] define the inverse rotation...
    if camera_id != 2
        new_sol[(1:2).+(9+(camera_id-3)*2)] = -sol[(1:2).+(9+(camera_id-3)*2)]
    else
        new_sol[10:13] = -sol[10:13]
    end
    return new_sol
end

function xyFlipSymmetry(sol)
    return vcat(-sol[1:6], sol[7:13])
end


sol = sols[rand(1:length(sols))]

# Verify Larsson symmetries
norm(F(larssonSymmetry(sol,2),params))
norm(F(larssonSymmetry(sol,3),params))
norm(F(larssonSymmetry(sol,4),params))

# Verify xyFlip symmetry
norm(F(xyFlipSymmetry(sol),params))

rqt1 = evaluateRQT(sol)
rqt2 = evaluateRQT(sols[rand(1:length(sols))])
M = [rqt1 rqt2]
svdvals(M)
