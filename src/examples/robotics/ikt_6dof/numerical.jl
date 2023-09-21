include("ikt_6dof.jl")
using HomotopyContinuation

function pose2vt(pose::Dict{String, Vector{Float64}})::Vector{Float64}
    return vcat(q2c(pose["q"]), pose["t"])
end

function joints2cs(joints::Vector{Float64}, offsets::Vector{Float64})::Vector{Float64}
    return vcat([cos(offsets[i] + joints[i]) for i in 1:6], [sin(offsets[i] + joints[i]) for i in 1:6])
end

function cs2joints(cs::Vector{Float64}, offsets::Vector{Float64})::Vector{Float64}
    return vcat([atan(cs[6+i], cs[i])-offsets[i] for i in 1:6])
end

function csSols2jointsSols(cs_sols::Vector{Vector{ComplexF64}}, offsets::Vector{Float64})::Vector{Vector{Float64}}
    return [cs2joints(real(cs), offsets) for cs in cs_sols if all(imag(cs) .< eps())]
end

function ikt_eqs_numerical(mechanism::Dict{String, Vector{Float64}})
    @var c[1:6], s[1:6], v[1:3], t[1:3]
    pose = [c2R(v) t; 0 0 0 1]

    # eqs = transformation_sym(mechanism, 5, 2, c, s) - inv_transformation_sym(mechanism, 2, 0, c, s) * pose * inv_transformation_sym(mechanism, 6, 5, c, s)
    # eqs = transformation_sym(mechanism, 6, 0, c, s) - pose

    T = transformation_sym(mechanism, 6, 0, c, s)
    R = T[1:3,1:3]
    trR = R[1,1]+R[2,2]+R[3,3]
    tran = T[1:3,4]

    eqs1 = xx2v(R'-R)-(1+trR)*v
    eqs2 = tran - t

    # eqs = vcat(M2V(eqs)...)
    eqs = vcat(eqs1, eqs2)
    for i in 1:6
        push!(eqs, c[i]^2 + s[i]^2 - 1)
    end
    return System(eqs; variables=vcat(c,s), parameters=vcat(v, t))
end

function incidence_6dof()
    @var c[1:6], s[1:6], v[1:3], t[1:3], a[1:6], d[1:6], α[1:6]
    
end

# Defining the mechanism...
mechanism = random_mechanism()
mechanism = Dict{String, Vector{Float64}}("theta offset" => randn(6),
                                          "d" => [350, 50, 0, 425, 0, 100],
                                          "a" => [50, 425, 0, 0, 0, 0],
                                          "alpha" => [-pi/2, 0, -pi/2, pi/2, -pi/2, 0])


# Setting the joints...
joints = randn(6)
cs = joints2cs(joints, mechanism["theta offset"])

# Computing the pose...
pose = fkt(mechanism, joints)
vt = pose2vt(pose)

# Forming the equations...
F = ikt_eqs_numerical(mechanism)

using DecomposingPolynomialSystems
grading = scaling_symmetries(F)

grading.U[1]

diff = differentiate(F.expressions, F.variables)
jac = exprDet(diff)

params = F.parameters
v1,v2,v3 = params[1],params[2],params[3]
t1,t2,t3 = params[4],params[5],params[6]

s = subs(jac, [v1,v2,v3,t1,t2,t3]=>[0,0,0,0,0,0])
println(s)

# Testing monodromy...
norm(F(cs, vt))
MR = monodromy_solve(F, [cs], vt, permutations=true)
cs_sols = solutions(MR)
joints_sols = csSols2jointsSols(cs_sols, mechanism["theta offset"])

# Galois group is S16
monodromy_permutations = filter_permutations_with_zeros(HomotopyContinuation.permutations(MR))
G = permutations_to_group(monodromy_permutations)
GAP.Globals.StructureDescription(G)

all_block_partitions(G)

# Test tracking...
joints2 = randn(6)
cs2 = joints2cs(joints2, mechanism["theta offset"])
pose2 = fkt(mechanism, joints2)
vt2 = pose2vt(pose2)

res = solve(F, cs, start_parameters=vt, target_parameters=vt2)
sols = solutions(res)[1]
min(svdvals([sols cs2])...)
