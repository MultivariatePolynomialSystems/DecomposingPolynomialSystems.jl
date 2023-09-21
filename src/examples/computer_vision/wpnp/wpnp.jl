gpu = false
if gpu
    include("../../compute_deck_transformations_gpu.jl")
else
    include("../../compute_deck_transformations_cpu.jl")
end

function wp_rot(q)
    return [q[1]^2 + q[2]^2 - q[3]^2 - q[4]^2      2*(q[2]*q[3] - q[1]*q[4])               2*(q[1]*q[3]+q[2]*q[4]);
            2*(q[1]*q[4] + q[2]*q[3])              q[1]^2 - q[2]^2 + q[3]^2 - q[4]^2       2*(q[3]*q[4] - q[1]*q[2])]
end

function squared_norm(v::Vector{Expression})
    return tr(transpose(v)*v)
end

function synthetize_solution()
    q = rand(ComplexF64, 4)
    R = wp_rot(q)
    a = rand(ComplexF64, 3)
    b = R*diagm(a)
    return ([q], vcat(reshape(a, 3), reshape(b, 6)))
end


@var q[1:4] a[1:3] b[1:2,1:3]
R = wp_rot(q)
M = reshape(R*diagm(a) - b, 6)
f = squared_norm(M)
eqs = reshape(differentiate([f], q), 4)

F = System(eqs; parameters=vcat(a, reshape(b, 6)), variables=q)
(x0, p0) = synthetize_solution()
norm(F(x0[1], p0))
deck_transformations = compute_deck_transformations(F, x0, p0, 3, tol=1e-5, expected_n_sols=32, param_dep=true)
