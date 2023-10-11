using DecomposingPolynomialSystems

function multidegreesFromTotalDegree(mds::Vector{Int}, n::Int, d::Int)::Vector{Vector{Int}}
    if n == 1
        return [vcat(mds, d)]
    end
    allMDs = Vector{Vector{Int}}([])
    for i in 0:d
        append!(allMDs, multidegreesFromTotalDegree(vcat(mds, [i]), n-1, d-i))
    end
    return allMDs
end

multidegreesFromTotalDegree(Vector{Int}([]), 4, 0)

function multidegreesUpToTotalDegree(n::Int, d::Int)::Vector{Vector{Int}}
    allMDs = Vector{Vector{Int}}([])
    for i in 0:d
        append!(allMDs, multidegreesFromTotalDegree(Vector{Int}([]), n, i))
    end
    return allMDs
end

multidegreesUpToTotalDegree(4, 2)

function remainderOfMultidegree(MD::Vector{Int}, A::Matrix{Int}, S::Vector{Int})
    @assert size(A, 1) == length(S)
    @assert length(MD) == size(A, 2)
    p = A*MD
    for i in 1:length(S)
        if S[i] > 0
            p[i] = mod(p[i], S[i])
        end
    end
    return p
end

MD = [1, 2, 3, 4]
A = [1 1 0 0; 0 0 1 1; 1 1 1 0]
S = [0, 0, 2]
remainderOfMultidegree(MD, A, S)

# rows of A describe the exponents of scaling symmetries
# elements of S whether it is finite (s[i] > 0) or infinite (s[i] == 0)
function partitionMultidegreesDynamic(MDs::Vector{Vector{Int}}, A::Matrix{Int}, S::Vector{Int})
    classes = Vector{Vector{Int}}([])
    rems = VV2M([remainderOfMultidegree(md, A, S) for md in MDs])
    rems .-= minimum(rems, dims=2) .- 1
    println(max(maximum(rems, dims=2)...))
    classesMatrix = zeros(Int, maximum(rems, dims=2)...)
    for i in 1:length(MDs)
        class_idx = classesMatrix[rems[:, i]...]
        if class_idx != 0
            push!(classes[class_idx], i)
        else
            classesMatrix[rems[:, i]...] = length(classes) + 1
            push!(classes, [i])
        end
        if mod(i, 100000) == 0
            println("MD id = ", i, " verified")
        end
    end
    return classes
end

function partitionMultidegrees(MDs::Vector{Vector{Int}}, A::Matrix{Int}, S::Vector{Int})
    classes = Vector{Vector{Int}}([])
    rems = [remainderOfMultidegree(md, A, S) for md in MDs]
    for i in 1:length(MDs)
        found = false
        for j in 1:length(classes)
            if rems[i] == rems[classes[j][1]]
                found = true
                push!(classes[j], i)
            end
        end
        if !found
            push!(classes, [i])
        end
        if mod(i, 100000) == 0
            println("MD id = ", i, " verified")
        end
    end
    return classes
end

MDs = multidegreesUpToTotalDegree(8, 2)
A = [1 1 0 0 1 1 0 0;
     0 0 1 1 0 0 1 1;
     1 1 1 1 0 0 0 0;
     0 0 0 0 1 1 1 1;
     1 1 0 0 0 0 0 0];
S = [0, 0, 0, 0, 2];
classes = partitionMultidegrees(MDs, A, S)

MDs = multidegreesUpToTotalDegree(8, 2)
A = [1 1 0 0 1 1 0 0;
     0 0 1 1 0 0 1 1;
     1 1 1 1 0 0 0 0;
     0 0 0 0 1 1 1 1];
S = [0, 0, 0, 0];
classes = partitionMultidegrees(MDs, A, S)

using HomotopyContinuation, LinearAlgebra
@var P[1:3,1:4,1:3] H[1:4,1:4] k[1:5] R[1:3,1:3,1:3] t[1:3,1:3] α[1:3]
Ps = [P[1:3,1:4,i] for i in 1:3]
Rs = [R[1:3,1:3,i] for i in 1:3]
ts = [t[1:3,i] for i in 1:3]
K = [k[1] k[2] k[3]; 0 k[4] k[5]; 0 0 1]
eqsP = vcat([M2V(Ps[i]-α[i]*K*[Rs[i] ts[i]]*H) for i in 1:3]...)
eqsRot = vcat(vcat([M2V(Rs[i]'*Rs[i]-eye(3)) for i in 1:3]...), [det(Rs[i])-1 for i in 1:3])
eqs = vcat(eqsP, eqsRot)
F = System(eqs; variables=vcat(M2V(P), M2V(H), k, M2V(R), M2V(t), α))
A, S = scalingSymmetries(F)
A = A[:, 1:length(M2V(P))]

# Monomial bases of homogeneous components
A = A[[5,6,7,8,9,10,12,13],:]
S = S[[5,6,7,8,9,10,12,13]]

MDs = multidegreesUpToTotalDegree(36, 6)
classes = partitionMultidegrees(MDs, A, S)
max(length.(classes)...)

a = [1.0, 2, 3]
b = [0, 2, 1]

@var x[1:3]
prod(x.^b)

function getSamples()::Vector{Vector{Float64}}
    Cs = [randn(3) for _ in 1:3]
    Rs = [c2R(c) for c in Cs]
    ts = [randn(3) for _ in 1:3]
    K = [randn() randn() randn(); 0 randn() randn(); 0 0 1]
    H = randn(4, 4)
    α = randn(3)
    return [VM2V([α[i]*K*[Rs[i] ts[i]]*H for i in 1:3])]
end

getSamples()

function interpolate(MDs::Vector{Vector{Int}}, MDclasses::Vector{Vector{Int}}, vars::Vector{Variable}, getSamples::Function)
    samples = getSamples()
    max_n_mons = max(length.(MDclasses)...)
    for i in 1:(div(max_n_mons, length(samples))+1)
        append!(samples, getSamples())
    end
    for MDclass in MDclasses
        if length(MDclass) == 1
            continue
        end
        n_samples = length(samples)
        n_mons = length(MDclass)
        A = zeros(typeof(samples[1][1]), length(samples), length(MDclass))
        monsMD = MDs[MDclass]
        for i in 1:n_samples
            A[i, :] = prod(samples[i].^VV2M(monsMD), dims=1)
        end
        N = nullspace(A)
        if size(N, 2) > 0
            println("Got Something!")
        end
    end
end

b = [1, 2, 3]
e = VV2M([[0, 0, 0], [1, 1, 1]])
b.^e
prod(b.^e, dims=1)

MDs = multidegreesUpToTotalDegree(36, 5)
classes = partitionMultidegrees(MDs, A, S)
interpolate(MDs, classes, VM2V(Ps), getSamples)

using DecomposingPolynomialSystems, AbstractAlgebra
@var T[1:2,1:2,1:2,1:2] P[1:4,1:4]
E = eye(4)
Ps = [[E[i,:] P[i,:]]' for i in 1:4]
eqs = [T[i1,i2,i3,i4] - (-1)^(i1+i2+i3+i4)*exprDet([transpose(P[i1,:,1]); transpose(P[i2,:,2]); transpose(P[i3,:,3]); transpose(P[i4,:,4])], expnd=false) for (i1,i2,i3,i4) in vec(reshape(collect(Iterators.product(ntuple(_ -> 1:2, 4)...)), 2^4, 1))]
F = System(eqs; variables=vcat(reshape(T, 16), reshape(P, 32)))
A, S = scalingSymmetries(F)
A = A[:, 1:16]

A = A[[1,2,3,6,10], :]
# hnf(matrix(ZZ, A))

# A = hnf(matrix(ZZ, A))
# A = Matrix{Int}(Matrix(A[1:5,:]))
S = S[[1,2,3,6,10]]

MDs = multidegreesUpToTotalDegree(16, 12)
classes = partitionMultidegrees(MDs, A, S)
max(length.(classes)...)

using DecomposingPolynomialSystems, LinearAlgebra

function eqs_5pp()
    @var R[1:3,1:3], t[1:3], α[1:5], β[1:5], x[1:3,1:5], y[1:3,1:5]
    eqs = vcat(reshape(transpose(R)*R - I, 9), [det(R) - 1])
    for i in 1:5
        append!(eqs, β[i]*y[:,i] - R*α[i]*x[:,i] - t)
    end
    return System(eqs; variables = vcat(reshape(R, 9), t, α, β), parameters = vcat(reshape(x, 15), reshape(y, 15)))
end

F = eqs_5pp()
A, S = scalingSymmetries(F)
A
S

MDs = multidegreesUpToTotalDegree(52, 4)
classes = partitionMultidegrees(MDs, A, S)
max(length.(classes)...)