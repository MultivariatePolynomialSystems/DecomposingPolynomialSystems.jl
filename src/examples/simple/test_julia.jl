function modify_array(M::Matrix{Int})
    M[1,1] = 100
end

M = [1 1 1; 2 2 2; 3 3 3]
modify_array(M)

M

function interpolate_symmetry_function(a::Vector{Int}, 
                                       b::Vector{Int})
    println(a)
    println(b)
end

interpolate_symmetry_function([1,2,3], [4,5,6])

mutable struct Monomial
    vars::Vector{Int}
end

function test(m::Monomial)
    a = [1,2,3]
    m.vars = a
    return a
end

m = Monomial([4,5,6])
a = test(m)

m.vars
a[1] = 5
a
m.vars

function prodMatrix(A::Matrix{Float64})
    return prod(A[1, :])
end

function prodSubmatrix(A::Matrix{Float64})
    return prod(view(A, 1, :))
end

A = rand(10, 100000)

using BenchmarkTools

@btime prodMatrix(A);
@btime prodSubmatrix(A);