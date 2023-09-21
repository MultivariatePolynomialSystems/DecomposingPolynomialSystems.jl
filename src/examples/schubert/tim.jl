using DecomposingPolynomialSystems, LinearAlgebra

@var K[1:4,1:8,1:4]
@var X[1:4,1:4]
@var Y[1:16,1:7]

XI = vcat(X,I)
constraint_matrices = [hcat(K[i,:,:],XI) for i=1:4]
n_mats = [C[1:7,1:7] for C in constraint_matrices]
s_mats = [C[2:8,2:8] for C in constraint_matrices]
e_mats = [C[2:8,1:7] for C in constraint_matrices]
w_mats = [C[1:7,2:8] for C in constraint_matrices]
Ms = vcat(n_mats,s_mats,w_mats,e_mats)
F = System(vcat([M2V(M*Y[i,:]) for (i, M) in enumerate(Ms)]...); parameters = K[:])

gr = scaling_symmetries(F, vcat(X[:], K[:]))
gr[1][2]

MDs = multidegrees_up_to_total_degree(144, 3)
classes = partition_multidegrees(MDs, gr)
max(length.(collect(values(classes)))...)

monodromy_solve(F)

xp0 = HomotopyContinuation.find_start_pair(HomotopyContinuation.fixed(F, compile=HomotopyContinuation.COMPILE_DEFAULT[]))
typeof(xp0)
F = run_monodromy(F)

@var a b c d x[1:2]
A = [a c; b d]
F = System(M2V(A*x))

gr = scaling_symmetries(F)
gr[1][2]