# a4 + b5 = log(1)
# a4 + b0 = log(1)
# a3 + b5 = log(2)
# a3 + b4 = log(4)
# a3 + b0 = log(2)
# a0 + b5 = log(3)
# a0 + b4 = log(6)
A = [0 0 1 0 0 1;
     0 0 1 1 0 0;
     0 1 0 0 0 1;
     0 1 0 0 1 0;
     0 1 0 1 0 0;
     1 0 0 0 0 1;
     1 0 0 0 1 0;
     0 0 0 0 0 1]
b = [log(1), log(1), log(2), log(4), log(2), log(3), log(6), log(1)]
using DecomposingPolynomialSystems
M = [A b]
M = rref(M)
S = M[1:6, :]
x = exp.(S[:,7])

A = [0 1 1;
     1 0 1;
     0 0 1]
b = [log(1), log(1), log(1)]
M = [A b]
M = rref(M)
x = exp.(M[:,4])

using DecomposingPolynomialSystems
using HomotopyContinuation
using LinearAlgebra

@var x y a b c
F = System([(x-a)^2*(y)^3 + a*(x-a)*(y)^2 + b*(x-a)*(y) + (x-a) + c, c*(x-a)^3*(y)^2 + (x-a)^3*(y) + b*(x-a)^2*(y) + a*(x-a)*(y) + 1]; variables=[x,y], parameters=[a,b,c])
F = run_monodromy(F, target_solutions_count=6, max_loops_no_progress=1000)
sample_system!(F, 5)


invs = compute_invariants(F, degree=1)
display(invs[1])

F.block_partitions[1]

A = zeros(ComplexF64, 15, 9)
for i in 1:5
     for j in 1:3
          M = ones(ComplexF64, 2, 6)
          M[:,1:2] = transpose(F.solutions[:, F.block_partitions[1][j], i])
          M[1,3:5] = F.parameters[:, i]
          M[2,3:5] = F.parameters[:, i]
          for k in 2:6
               A[(i-1)*3+j,k-1] = M[1,1]*M[2,k]-M[2,1]*M[1,k]
          end
          for k in 3:6
               A[(i-1)*3+j,6+k-3] = M[1,2]*M[2,k]-M[2,2]*M[1,k]
          end
     end
end
N = sparsify(nullspace(A), 1e-5)