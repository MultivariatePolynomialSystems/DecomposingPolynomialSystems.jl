using LinearAlgebra, BenchmarkTools

function eig20(n::Int)
    for i in 1:n
        e = eigvals(randn(20, 20))
    end
end

function eig10(n::Int)
    for i in 1:n
        e = eigvals(randn(10, 10))
        for j in 1:10
            h = eigvals(randn(2, 2))
        end
    end
end

n = 1000
@btime eig20(n)
@btime eig10(n)