using HomotopyContinuation

n = 1
d = 3

function monomials_exponents(n, d; affine::Bool)
    if affine
        pred = x -> sum(x) ≤ d
    else
        pred = x -> sum(x) == d
    end
    E = map(Iterators.filter(pred, Iterators.product(Iterators.repeated(0:d, n)...))) do e
        collect(e)
    end

    sort!(E, lt = td_order)
    E
end

function td_order(x, y)
    sx = sum(x)
    sy = sum(y)
    sx == sy ? x > y : sx > sy
end

monomials_exponents(n, d, affine=true)

collect(Iterators.product(Iterators.repeated(0:d, n)...))

# In the Iterators.product(Iterators.repeated(0:d, n)...)
# part the code will generate (d+1)^n n-tuples and among them pick the ones that fullfil
# the total degree condition (≤ d or == d). The problem is that with large enough n (or d)
# for example n = 20, d = 3, the part E = ... becomes a bottleneck since 4^20 is a huge number.

# I suggest the reccursion solution to this problem.

function multidegrees_from_total_degree(md::Vector{Int}, n::Int, d::Int)::Vector{Vector{Int}}
    if n == 1
        return [vcat(md, d)]
    end
    mds = Vector{Vector{Int}}([])
    for i in 0:d
        append!(mds, multidegrees_from_total_degree(vcat(md, i), n-1, d-i))
    end
    return mds
end

function multidegrees_up_to_total_degree(n, d)::Vector{Vector{Int}}
    if n ≤ 0
        return []
    end
    mds = Vector{Vector{Int}}([])
    for i in 0:d
        append!(mds, multidegrees_from_total_degree(Vector{Int}([]), n, i))
    end
    return mds
end

multidegrees_up_to_total_degree(n, d)