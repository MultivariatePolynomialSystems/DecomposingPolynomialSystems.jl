using DecomposingPolynomialSystems, AbstractAlgebra

function divides(a::Int, b::Int)::Bool
    return b == a*div(b,a)
end

function allTuplesGivenWeightedSum(elems::Vector{Int}, A::Vector{Int}, b::Int)::Vector{Vector{Int}}
    if length(A) == 1
        if divides(A[1], b)
            return [vcat(elems, div(b, A[1]))]
        else
            return []
        end
    end
    allTuples = Vector{Vector{Int}}([])
    for i in 0:div(b, A[1])
        append!(allTuples, allTuplesGivenWeightedSum(vcat(elems, [i]), A[2:end], b-i*A[1]))
    end
    return allTuples
end

allTuplesGivenWeightedSum(Vector{Int}([]), [1,1,1,1,1], -4)

# For a linear diophantic equation to have finite set of solutions
# in non-negative integers a necessary condition is A > 0 and b >= 0
function solveNonNegative(A::Vector{Int}, b::Int)::Vector{Vector{Int}}
    if (!all(A .> 0) || b < 0 || !divides(gcd(A), b))
        return []
    end
    return allTuplesGivenWeightedSum(Vector{Int}([]), A, b)
end

A = [3,1,2]
b = 3
solveNonNegative(A, b)

function solveNonNegativeTriangular(A::Matrix{Int}, b::Vector{Int}, sol::Vector{Int})::Vector{Vector{Int}}
    if size(A, 1) == 1
        sols = solveNonNegative(M2V(A), b[1])
        return [vcat(s, sol) for s in sols]
    end
    A₁, b₁ = A[end,findfirst(a->a!=0, A[end,:]):end], b[end]
    n_vars = length(A₁)
    partialSols = solveNonNegative(A₁, b₁)
    allSols = Vector{Vector{Int}}([])
    for partialSol in partialSols
        A₂ = A[1:end-1,1:end-n_vars]
        b₂ = b[1:end-1]-A[1:end-1,(end-n_vars+1):end]*partialSol
        append!(allSols, solveNonNegativeTriangular(A₂, b₂, vcat(partialSol, sol)))
    end
    return allSols
end

A = [1 1 0 1; 
     0 0 1 1]
d = 1
b = d*[1; 2]
solveNonNegativeTriangular(A, b, Vector{Int}([]))

function toJM(M)
    return Matrix{Int}(Matrix(M))
end

# The rows of A encode the exponents of infinite scaling symmetries,
# d is the upper bound for the degree of monomials
function monomialBasesOfHomogeneousComponents(A::Matrix{Int}, d::Int)
    n_vars = size(A, 2)
    A = [A; ones(Int, 1, size(A, 2))]
    A = matrix(ZZ, [A -1*eye(size(A, 1))])
    H = hnf(A)
    display(H)
    A = H[:,1:end-1]
    c = -d*H[:,end]
    j = 0
    for i in size(A,1):-1:1
        if A[i,1:n_vars] != matrix(ZZ, zeros(Int, 1, n_vars))
            j = i
            break
        end
    end
    A, M = toJM(A[1:j,:]), toJM(A[j+1:end,n_vars+1:end])
    b, m = M2V(toJM(c[1:j,1])), M2V(toJM(c[j+1:end,1]))
    solsM = solveNonNegativeTriangular(M, m, Vector{Int}([]))
    nM = size(M, 2)
    monBases = Vector{Vector{Vector{Int}}}([])
    for solM in solsM
        Aₛ = A[:,1:end-nM]
        bₛ = b-A[:,(end-nM+1):end]*solM
        expSols = solveNonNegativeTriangular(Aₛ, bₛ, Vector{Int}([]))
        push!(monBases, expSols)
    end
    return monBases
end

A = [1 1 0 0 1 1 0 0;
     0 0 1 1 0 0 1 1;
     1 1 1 1 0 0 0 0;
     0 0 0 0 1 1 1 1]
d = 2
# monBases = monomialBasesOfHomogeneousComponents(A, d)

using DecomposingPolynomialSystems, HomotopyContinuation, AbstractAlgebra

function M2ColDiffs(M::Matrix{Int})
    M = M - M[:,1]*ones(Int, 1, size(M,2))
    return M[:,2:end]
end

function modV(V::Vector{BigInt}, n::BigInt)
    return [mod(v, n) for v in V]
end

function scalingSymmetries(F::System)
    Es = [Matrix{Int}(exponents_coefficients(f, variables(F))[1]) for f in F.expressions]
    Es = [M2ColDiffs(E) for E in Es]
    A = matrix(ZZ, VM2M(Es))
    @assert size(A, 1) <= size(A, 2)
    S, T, _ = snf_with_transform(A)
    S, T = Matrix(S), Matrix(T)
    i = 1
    while i <= size(S, 1) && S[i,i] == 1
        i += 1
    end
    j = i
    while j <= size(S,1) && S[j,j] != 0
        T[j, :] = modV(T[j,:], S[j,j])
        j += 1
    end
    return Matrix{Int}(Matrix(T[i:end, :])), [Int(S[j,j]) for j in i:size(S, 1)]
end

@var x a b
F = System([x^2+a*b+b, x^2 + a])
T, S = scalingSymmetries(F)
T
S

@var P[1:3,1:4,1:3] H[1:4,1:4] k[1:5] R[1:3,1:3,1:3] t[1:3,1:3] α[1:3]
Ps = [P[1:3,1:4,i] for i in 1:3]
Rs = [R[1:3,1:3,i] for i in 1:3]
ts = [t[1:3,i] for i in 1:3]
K = [k[1] k[2] k[3]; 0 k[4] k[5]; 0 0 1]
eqsP = vcat([M2V(Ps[i]-α[i]*K*[Rs[i] ts[i]]*H) for i in 1:3]...)
eqsRot = vcat(vcat([M2V(Rs[i]'*Rs[i]-eye(3)) for i in 1:3]...), [det(Rs[i])-1 for i in 1:3])
eqs = vcat(eqsP, eqsRot)
F = System(eqs; variables=vcat(VM2V(Ps), M2V(H), k, VM2V(Rs), M2V(t), α))
T, S = scalingSymmetries(F)
T = T[:, 1:length(VM2V(Ps))]

# Monomial bases of homogeneous components
A = T[[5,6,7,8,9,10,12,13],:]
hnf(matrix(ZZ, A))

monomialBasesOfHomogeneousComponents(A, 2)