# Using symmetries in the eigenvalue method for polynomial systems
# Corless, Gattermann, Kotsireas
using DecomposingPolynomialSystems
using HomotopyContinuation
using LinearAlgebra

@var x a b
F = System([a*x^4 + b*x^3 + x^2 + b*x + a]; variables=[x], parameters=[a,b])
F = run_monodromy(F)
f = x+1/x # invariant under the symmetry x --> 1/x
# B = [x^2-1/x^2, x-1/x, x+1/x, 1] # bases of isotypic components
B = [x^3, x^2, x+1/x, 1]
is_basis(F, B)
A = multiplication_matrix(F, f, B, 1) # multiplication matrix w.r.t. symmetry-adapted basis

A = A[3:4,3:4]
E = eigen(A)
V = p2a(E.vectors) # evaluations of x+1/x at solutions
[x₀+1/x₀ for x₀ in F.solutions[1, :, 1]]

x+1/x = t
#________________________________________________________________________________#

using DecomposingPolynomialSystems
using HomotopyContinuation
using LinearAlgebra

@var x y a b c
F = System([(x-a)^2*(y)^3 + a*(x-a)*(y)^2 + b*(x-a)*(y) + (x-a) + c, c*(x-a)^3*(y)^2 + (x-a)^3*(y) + b*(x-a)^2*(y) + a*(x-a)*(y) + 1]; variables=[x,y], parameters=[a,b,c])
F = run_monodromy(F, target_solutions_count=6, max_loops_no_progress=1000)
sample_system!(F, 200)
sample_system_once!(F, ComplexF64.([-1, 1, -1]))

B = [x^2, x*y, y^2, x, y, 1]
f₁ = (x-a)/y
f₂ = y/(x-a)
express_in_basis(F, f₁, B, 4)

G = permutations_to_group(F.monodromy_permutations)
GAP.Globals.MinimalGeneratingSet(G)
group_structure(G)
centralizer(G)

symmetries = compute_symmetries!(F, degree=2)
B₁ = monomial_basis(F)

B = [x^2, x*y, y^2, x, y, 1]
instance_id, a₀, b₀, c₀ = 201, -1, 1, -1
M1 = multiplication_matrix(F, Expression(x), B, instance_id) # mult by x in B₁
M2 = multiplication_matrix(F, Expression(y), B, instance_id) # mult by y in B₂

M3 = sparsify(M1 + inv(M2)*(eye(6) - M2), 1e-5)

M4 = sparsify(inv(M2)*(M1 + eye(6)), 1e-5)
M5 = sparsify(inv(M1 + eye(6))*M2, 1e-5)

T = zeros(ComplexF64, 6, 6)
T[:, 1] = [1, 0, 0, 0, 0, 0]
T[:, 2] = [0, 1, 0, 0, 0, 0]
T[:, 3] = [0, 0, 1, 0, 0, 0]
T[:, 4] = M4[:, 6]
T[:, 5] = M5[:, 6]
T[:, 6] = [0, 0, 0, 0, 0, 1]
T

M6 = sparsify(7*inv(T)*M3*T, 1e-5)

M₃ = sparsify(M₁ + (eye(6)+a₀*M₂)inv(M₂), 1e-5) # mult by x + (1+ay)/y in B₁
multiplication_matrix(F, x+(1+a₀*y)/y, B₁, instance_id)

M₄ = sparsify((M₁-a₀*eye(6))*inv(M₂), 1e-5) # mult by (x-a)/y in B₁
multiplication_matrix(F, (x-a₀)/y, B₁, instance_id)

M₅ = sparsify(M₂*inv(M₁-a₀*eye(6)), 1e-5) # mult by y/(x-a) in B₁
multiplication_matrix(F, y/(x-a₀), B₁, instance_id)

T = zeros(ComplexF64, 6, 6)
T[:, 1] = [1, 0, 0, 0, 0, 0]
T[:, 2] = M₄[:, 1]
T[:, 3] = M₅[:, 1]
T[:, 4] = [0, 1, 0, 0, 0, 0]
T[:, 5] = [0, 0, 1, 0, 0, 0]
T[:, 6] = [0, 0, 0, 0, 1, 0]
T

B₂ = [1, x/y, y/x, x, y, x*y]
M₆ = sparsify(inv(T)*M₃*T, 1e-5) # mult by x + (1+ay)/y in B₂

B = [(x-a₀)/y, y/(x-a₀), 1, x, y, x*y]
M = multiplication_matrix(F, x+(1+a₀*y)/y, B, instance_id)

A = transpose(M₆[1:3,1:3])
E = eigen(A)
λ = E.values
V = E.vectors
V = V./reshape(V[1,:], 1, 3)
V = V[2:3,:]

F₁ = System([x*y+(a₀-λ[1])*y+1, x-V[1,1]*y-a₀])
S₁ = solutions(solve(F₁))

F₂ = System([x*y+(a₀-λ[2])*y+1, x-V[1,2]*y-a₀])
S₂ = solutions(solve(F₂))

F₃ = System([x*y+(a₀-λ[3])*y+1, x-V[1,3]*y-a₀])
S₃ = solutions(solve(F₃))

F.solutions[:, :, 201]


using DecomposingPolynomialSystems

A = Matrix{Float64}([0  -1 0  -1 -2 0  0  1  -2 -1 1  -1 0  -1 0  0;
0  0  -1 0  -1 -2 -1 0  1  -2 0  1  -1 0  -1 0;
0  0  1  0  1  1  0  0  1  0  -1 2  0  1  0  0;
0  -1 1  0  -1 1  0  0  -2 1  -1 -1 0  0  0  0;
0  0  0  0  -1 0  0  0  -1 -2 1  -2 -1 1  -1 -1;
0  0  0  0  0  1  1  0  0  1  0  1  0  -1 2  1;
0  1  0  0  3  0  0  0  -1 3  2  -2 1  1  1  -1;
0  0  0  -1 0  0  0  0  -4 0  -1 -1 -2 0  0  0;
0  0  0  0  -1 1  1  0  0  -1 0  -2 1  -1 -1 0;
-1 0  0  -5 -1 0  0  -1 -8 -3 1  -7 -2 1  -2 -2;
1  0  0  5  0  2  4  2  7  -1 5  4  -1 0  5  1])
rref(A, 1e-5)

using LinearAlgebra
rank(A)

# a₂ --> a₂ + c₂ + d₂ - b₃ - b₄
# b₂ --> b₂
# c₂ --> b₃
# d₂ --> b₄
# a₃ --> a₃ + b₃ + d₃ - c₂ - c₄
# b₃ --> c₂
# c₃ --> c₃
# d₃ --> c₄
# a₄ --> a₄ + b₄ + c₄ - d₂ - d₃
# b₄ --> d₂
# c₄ --> d₃
# d₄ --> d₄

# c₂ = b₃
# d₂ = b₄
# d₃ = c₄

# Invariant monomials
# p₁₁