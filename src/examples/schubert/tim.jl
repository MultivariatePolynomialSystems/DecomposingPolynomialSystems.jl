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
    return System(dets; parameters = K[:])
end

F₁ = F1()
scalings = scaling_symmetries(F₁, vcat(X[:], K[:]))
@show scalings.action[1]

MDs = multidegrees_up_to_total_degree(144, 3)
classes = partition_multidegrees(MDs, scalings.grading)
max(length.(collect(values(classes)))...)

F₂ = F2()
monodromy_solve(F₂)
monodromy_solve(F₂; permutations=true)
F₂ = run_monodromy(F₂)
