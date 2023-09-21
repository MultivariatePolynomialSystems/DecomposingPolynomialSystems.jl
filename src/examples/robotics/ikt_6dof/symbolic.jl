include("ikt_6dof.jl")
using Singular, GroebnerBasis
using LinearAlgebra
import HomotopyContinuation: @var, differentiate

function ikt_eqs_symbolic(mechanism::Dict{String, Float64}, pose::Dict{String, Vector{Float64}}, tol::Float64)
    @var c[1:6], s[1:6]
    str_vars = [string(var) for var in vcat(c, s)]
    R, vars = PolynomialRing(QQ, str_vars)
    c, s = vars[1:6], vars[7:12]

    rat_pose = Matrix{Any}(rational_pose(pose, tol))
    rat2nQ!(rat_pose)
    rat_pose = Matrix{n_Q}(rat_pose)

    rat_mechanism = rational_mechanism(mechanism, tol)
    eqs = transformation_sym(rat_mechanism, 5, 2, c, s) - inv_transformation_sym(rat_mechanism, 2, 0, c, s) * rat_pose * inv_transformation_sym(rat_mechanism, 6, 5, c, s)
    # eqs = transformation_sym(rat_mechanism, 6, 0, c, s) - rat_pose
    eqs = vcat(M2L(eqs)...)
    for i in 1:6
        push!(eqs, c[i]^2 + s[i]^2 - 1)
    end
    return Ideal(R, eqs)
end

offsets = [0, -pi/2, 0, 0, -pi/2, 0]
mechanism = Dict{String, Float64}("d1" => 0.5, "a1" => 0.15, "theta1 offset" => offsets[1], "alpha1" => -pi/2,
                                  "d2" => 0, "a2" => 0.614, "theta2 offset" => offsets[2], "alpha2" => 0,
                                  "d3" => 0, "a3" => 0.2, "theta3 offset" => offsets[3], "alpha3" => -pi/2,
                                  "d4" => 0.64, "a4" => 0, "theta4 offset" => offsets[4], "alpha4" => pi/2,
                                  "d5" => 0, "a5" => 0.03, "theta5 offset" => offsets[5], "alpha5" => -pi/2,
                                  "d6" => 0.2, "a6" => 0.0, "theta6 offset" => offsets[6], "alpha6" => 0)
joints = [pi, pi, 0, pi/2, pi/2, 0]
pose = fkt(mechanism, joints)
tol = 1e-5

eqs = ikt_eqs_symbolic(mechanism, pose, tol)
gb = f4(eqs)
