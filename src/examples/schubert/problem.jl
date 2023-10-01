using DecomposingPolynomialSystems, LinearAlgebra, HomotopyContinuation

@var K[1:4,1:8,1:4]
@var X[1:4,1:4]
@var Y[1:16,1:7]

function F1()
    XI = vcat(X,I)
    constraint_matrices = [hcat(K[i,:,:],XI) for i=1:4]
    n_mats = [C[1:7,1:7] for C in constraint_matrices]
    s_mats = [C[2:8,2:8] for C in constraint_matrices]
    e_mats = [C[2:8,1:7] for C in constraint_matrices]
    w_mats = [C[1:7,2:8] for C in constraint_matrices]
    Ms = vcat(n_mats,s_mats,w_mats,e_mats)
    return System(vcat([M2V(M*Y[i,:]) for (i, M) in enumerate(Ms)]...))
end

function F2()
    XI = vcat(X,I)
    constraint_matrices = [hcat(K[i,:,:],XI) for i=1:4]
    n_dets = [det(C[1:7,1:7]) for C in constraint_matrices]
    s_dets = [det(C[2:8,2:8]) for C in constraint_matrices]
    e_dets = [det(C[2:8,1:7]) for C in constraint_matrices]
    w_dets = [det(C[1:7,2:8]) for C in constraint_matrices]
    dets = vcat(n_dets,s_dets,w_dets,e_dets)
    return System(dets; variables=X[:], parameters = K[:])
end

function fabricate_sample()
    X0 = randn(CC, 4, 4)
    K0_parts = [randn(CC, 8, 2) for i=1:4]
    K0s = [hcat(K0, hcat(vcat(X0, I), K0) * randn(CC, 6, 2)) for K0 in K0_parts]
    M2V(X0), M2V([K0s[i][j, k] for i=1:4, j=1:8, k=1:4])
end

F₁ = F1()
scalings = scaling_symmetries(F₁, vcat(X[:]))

F₂ = F2()
xp0 = fabricate_sample()
F₂ = run_monodromy(F₂, xp0)

deck = symmetries_fixing_parameters_graded!(F₂, scalings; degree_bound=6)
deck

mds = multidegrees_up_to_total_degree(144, 3)
classes = partition_multidegrees(mds, scalings.grading)
max(length.(collect(values(classes)))...)
