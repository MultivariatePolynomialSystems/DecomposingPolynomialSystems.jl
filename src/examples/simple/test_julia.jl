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

struct Foo
    bar
    baz::Vector{Int}
end

foo = Foo("Hello", [5])
typeof(foo.bar)

push!(foo.baz, 1)

abstract type AlgebraicGroup end

struct DeckTransformationGroup <: AlgebraicGroup
    vars::Vector{Int}
end

DTG = DeckTransformationGroup([1, 2, 3])

function print_vars(G::AlgebraicGroup)
    @show G.vars
end

print_vars(DTG)

struct Polar{T<:Real} <: Number
    r::T
    Θ::T
end
Polar(r::Real,Θ::Real) = Polar(promote(r,Θ)...)
Base.show(io::IO, z::Polar) = print(io, z.r, " * exp(", z.Θ, "im)")
Base.show(io::IO, ::MIME"text/plain", z::Polar{T}) where{T} = print(io, "Polar{$T} complex number:\n   ", z)

Polar(3, 4.0)
[Polar(3, 4.0), Polar(4.0,5.3)]

Base.show(io::IO, ::MIME"text/html", z::Polar{T}) where {T} = println(io, "<code>Polar{$T}</code> complex number: ", z.r, " <i>e</i><sup>", z.Θ, " <i>i</i></sup>")

show(stdout, "text/html", Polar(3.0,4.0))

# ----------------------------------------------------------------------------
using BenchmarkTools
M2V(M::AbstractMatrix) = M[:]
M = randn(1000, 1000)
@btime v = M2V(M);
@btime v = view(M, :);

# ----------------------------------------------------------------------------
using BenchmarkTools, LinearAlgebra
method1(v, M; tol=1e-5) = findfirst(i->norm(M[:,i]-v)<tol, axes(M,2))
method2(v, M; tol=1e-5) = findfirst(i->norm(view(M,:,i)-v)<tol, axes(M,2))
M2VV(M::AbstractMatrix) = [M[:,i] for i in axes(M, 2)]
method3(v, M; tol=1e-5) = findfirst(x->norm(x-v)<tol, M2VV(M))
M = randn(1000, 1000)
v = M[:,end]
@btime method1(v, M);
@btime method2(v, M);
@btime method3(v, M);

# ----------------------------------------------------------------------------
M2V(M::AbstractMatrix) = M[:]
M = randn(1000, 1000)
using BenchmarkTools
@btime v = M2V(M);
@btime v = view(M, :);
v = view(M, :)

M = randn(1000, 1000)
@views v = M2V(M)
v .= 0
M

M = randn(1000, 1000)
v = view(M, :)
v .= 0
M

M = randn(1000, 1000)
@views v = M[:]
v .= 0
M

M = randn(100, 100)
v = view(M, :)
a = copy(v)
v .= 0
M
a

M2VV(
    M::AbstractMatrix;
    dim::Integer,
    cp::Bool
) = [ifelse(cp, copy(selectdim(M, dim, i)), selectdim(M, dim, i)) for i in axes(M, dim)]
M = randn(2, 3)
v = M2VV(M; dim=2, cp=false)

# ----------------------------------------------------------------------------
using BenchmarkTools
function func1(M)
    return copy(view(M, :))
end
function func2(M)
    return M[:]
end
function func3(M)
    A = zeros(eltype(M), length(M))
    return copyto!(A, view(M, :))
end
function func4(M)
    return vec(M)
end
function func5(M)
    return view(M, :)
end
M = randn(1000, 1000)
@btime v = func1(M);
@btime v = func2(M);
@btime v = func3(M);
@btime v = func4(M);
@btime v = func5(M);

M[10000]

x = 0
while x < 10
    x += 1
end

# ----------------------------------------------------------------------------
Tuple{Vector{Int}, Vector{Int}} <:  Tuple{Vector{<:Integer}, Vector{<:Integer}}

Dict{Vector{Int}, Vector{Int}} <:  Dict{Vector{<:Integer}, Vector{<:Integer}}
Dict{Vector{Int}, Vector{Int}} <:  Dict{Vector{<:Integer}, Vector{Int}}
Dict{Vector{Int}, Vector{Int}} <:  Dict{Vector{Int}, Vector{<:Integer}}
Dict{Vector{Int}, Vector{Int}} <:  Dict{Vector{T}, Vector{T}} where T <: Integer
Dict{Vector{Int}, Vector{Int}} <:  Dict{Vector{T}, Vector{S}} where {T <: Integer, S <:Integer}

d = Dict{Vector{<:Integer}, Vector{<:Integer}}()
d[[1]] = [1]
typeof(d)

typeof(Tuple{Vector{T}, Vector{T}} where T <: Integer)
typeof(Tuple{Vector{<:Integer}, Vector{<:Integer}})
t = Tuple{Vector{<:Integer}, Vector{<:Integer}}(([1], [2]))

supertypes(Dict{Vector{<:Integer}, Vector{<:Integer}})
t = Dict{Vector{T}, Vector{T}} where T <: Integer
fieldnames(t)



typeof(Array{Array{T, 1} where T, 1})
typeof(Array{Array{T, 1}, 1} where T)

a = Array{Array{T, 1} where T, 1}()
push!(a, [1,2,3])
push!(a, [1.0])

convert(Vector{Int32}, [1,2,3])

a = ([1,2,3], [4,5,6])

struct Test
    n::Integer
end

t = Test(1)
typeof(t.n)

[i+j for i = 1:3 for j = 1:3]

function test(M::AbstractMatrix{T})::Matrix{T} where {T<:Number}
    A = zeros(T, 2, 2)
    A[1,1] = M[1,1]
    return A
end

M = ones(Complex, 3, 3)
test(M)