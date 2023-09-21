include("../compute_factorizing_maps.jl")
include("../../Gauss-Jordan.jl")
include("../../utils_cpu.jl")

@var x p


# F = System([x^4 + x^2 + p; parameters = [p]) # branch point at infinity
# q*x^4 + q*x^2*y^2 + p*y^4 = 0; q = 0 --> p = 1 -- > y^4 = 0 --> y = 0
# [p : q] = [1 : 0], [x : y] = { [1 : 0]x4 } --> gives 4-cyclic permutation
# [p : q] = [0 : 1], [x : y] = { [0 : 1]x2, [i : 1], [-i : 1] }
# [p : q] = [1 : 4] --> 4x^4 + 4x^2y^2 + y^4, y = 1 --> 4x^4 + 4x^2 + 1 = 0 --> D = 16 - 16 = 0,


F = System([x^4 + x^2 + p]; parameters = [p])
x0 = Vector{Vector{ComplexF64}}([[1]])
p0 = Vector{ComplexF64}([-2])
MR = monodromy_solve(F, x0, p0, permutations=true, max_loops_no_progress=50)
perms = HomotopyContinuation.permutations(MR)
G = permutations_to_group(perms)
GAP.Globals.StructureDescription(G)
GAP.Globals.Order(G)
centralizer(G)

all_block_partitions(G)



d = 2 # degree of the factorizing map
m = 2*(d+1) # number of unknown coefficients of the factorizing map
n_mons = d*(d+1) # number of binomials occuring in the equations
n = Int(n_mons/2) # number of unknowns of a linear system to solve (since 2i-th columns are -1 multiples of the previous ones)
X = zeros(ComplexF64, n, 2)
X[:, 1] = rand(ComplexF64, size(X, 1))
X[:, 2] = -X[:, 1]

M = zeros(ComplexF64, n, n)
for i = 1:n
    x = X[i, 1]
    y = X[i, 2]
    M[i, :] = [x^2*y-x*y^2, x^2-y^2, x-y]
end
M

tol = 1e-5
coeffs = rref(Matrix{ComplexF64}(transpose(nullspace(M))), tol)
coeffs = sparsify(coeffs, tol)


comesFromConstantMap(coeffs[1, :])
comesFromConstantMap(coeffs[2, :])
comesFromConstantMap(coeffs[3, :])
comesFromConstantMap(coeffs[4, :])

using Singular
R, vars = PolynomialRing(QQ, ["a0", "a1", "a2", "b0", "b1", "b2", "w1", "w2", "w3", "w4", "w5", "w6"])
a = vars[1:3]
b = vars[4:6]
w = vars[7:12]
gens = Vector{spoly{n_Q}}([])
k = 1
for i = 1:length(a)
    for j = 1:length(b)
        if i != j
            push!(gens, w[k]-a[i]*b[j])
            k += 1
        end
    end
end
I = Ideal(R, gens)
J = eliminate(I, vcat(a, b)...)

R, vars = PolynomialRing(QQ, ["a0", "a1", "a2", "a3", "b0", "b1", "b2", "b3", "w1", "w2", "w3", "w4", "w5", "w6", "w7", "w8", "w9", "w10", "w11", "w12"])
a = vars[1:4]
b = vars[5:8]
w = vars[9:20]
gens = Vector{spoly{n_Q}}([])
k = 1
for i = 1:length(a)
    for j = 1:length(b)
        if i != j
            push!(gens, w[k]-a[i]*b[j])
            k += 1
        end
    end
end
I = Ideal(R, gens)
J = eliminate(I, vcat(a, b)...)

using SmithNormalForm, LinearAlgebra
n = 3 # number of monomials
A = zeros(Int64, 2*n, n*(n-1)+1)
A[1, 1] = 1
k = 2
for i = 1:n
    for j = 1:n
        if i != j
            A[i, k] = 1
            A[n + j, k] = 1
            k += 1
        end
    end
end
A
F = smith(A)
diagm(F)


# Checking for algebraic independence
using HomotopyContinuation
@var x,y
f = x+y
g = (x+y)^2 + x + y
J = differentiate([f,g],[x,y])
rref(J)

# Elimination ideal
using Singular
using HomotopyContinuation
R, vars = PolynomialRing(QQ, [""])

@var y[1:3]
