using DecomposingPolynomialSystems, HomotopyContinuation, LinearAlgebra

function fabricateSolution()
    cs = [zeros(3), randn(ComplexF64, 3), randn(ComplexF64, 3)]
    Rs = [c2R(c) for c in cs]
    ts = [zeros(3), randn(ComplexF64, 3), randn(ComplexF64, 3)]
    Ps = [[Rs[i] ts[i]] for i in 1:3]
    X = randn(ComplexF64, 4, 3)
    X = p2a([X X[:,2:3]*randn(ComplexF64, 2)])
    xs = [Ps[i]*a2p(X) for i in 1:3]
    αs = [reshape(xs[i][3,:], 1, 4) for i in 1:3]
    xs = [xs[i]./αs[i] for i in 1:3]
    s = αs[1][1]
    αs = [αs[i]/s for i in 1:3]
    as = [-M2V(p2a(nullspace(xs[i][:,2:4]))) for i in 1:3]
    # println([as[i][1]+as[i][2] for i in 1:3])
    as = [as[i][1] for i in 1:3]
    # println([xs[i][:,2:3]*as[i]-xs[i][:,4] for i in 1:3])
    return (Vector{ComplexF64}(vcat(αs[1][2:4], reshape(αs[2], 4), reshape(αs[3], 4))), vcat(reshape(xs[1][1:2,1:3], 6), reshape(xs[2][1:2,1:3], 6), reshape(xs[3][1:2,1:3], 6), as))
end

fabricateSolution()

function fabricateSolution2()
    cs = [zeros(3), randn(ComplexF64, 3), randn(ComplexF64, 3)]
    Rs = [c2R(c) for c in cs]
    ts = [zeros(3), randn(ComplexF64, 3), randn(ComplexF64, 3)]
    Ps = [[Rs[i] ts[i]] for i in 1:3]
    X = randn(ComplexF64, 4, 3)
    a = randn(ComplexF64, 2)
    X = [X X[:,2:3]*a]
    x = [Ps[i]*X for i in 1:3]
    return (ones(ComplexF64, 11), vcat(M2V(x[1][:,1:3]), M2V(x[2][:,1:3]), M2V(x[3][:,1:3]), a, a, a))
end

xp0 = fabricateSolution2()

function eqs_3v4p()
    @var α₁[2:4], α₂[1:4], α₃[1:4] x[1:2,1:3,1:3], a[1:3]
    a₁ = [1; α₁]
    xa = [[x[:,:,i] x[:,2:3,i]*[a[i]; 1-a[i]]] for i in 1:3]
    eqs = []
    for i in 1:4
        for j in (i+1):4
            dA = a₁[i]*[xa[1][:,i]; 1] - a₁[j]*[xa[1][:,j]; 1]
            dB = α₂[i]*[xa[2][:,i]; 1] - α₂[j]*[xa[2][:,j]; 1]
            dC = α₃[i]*[xa[3][:,i]; 1] - α₃[j]*[xa[3][:,j]; 1]
            append!(eqs, transpose(dA)*dA - transpose(dB)*dB)
            append!(eqs, transpose(dA)*dA - transpose(dC)*dC)
        end
    end
    return System(eqs; variables = vcat(α₁, α₂, α₃), parameters = vcat(reshape(x, 18), a))
end

F = eqs_3v4p()
xp0 = find_start_pair(F)
norm(F(xp0[1], xp0[2]))
J = jac(F, xp0)
svdvals(J)

F = run_monodromy(F, xp0, max_loops_no_progress=1000)