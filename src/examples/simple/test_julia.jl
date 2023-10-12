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

findfirst(sol::Vector{CC}, sols::Matrix{CC}, tol::Float64) = findfirst(x->norm(x-sol)<tol, M2VV(sols))