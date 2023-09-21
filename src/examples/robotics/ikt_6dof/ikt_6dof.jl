include("../../../utils/utils.jl")
import Singular: spoly, n_Q

function q2r(q::Vector{<:Number})::Matrix{<:Number}
    return 1/(q[1]^2+q[2]^2+q[3]^2+q[4]^2)*[q[1]^2+q[2]^2-q[3]^2-q[4]^2 2*(q[2]*q[3]-q[1]*q[4]) 2*(q[2]*q[4]+q[1]*q[3]);
                     2*(q[2]*q[3]+q[1]*q[4]) q[1]^2-q[2]^2+q[3]^2-q[4]^2 2*(q[3]*q[4]-q[1]*q[2]);
                     2*(q[2]*q[4]-q[1]*q[3]) 2*(q[3]*q[4]+q[1]*q[2]) q[1]^2-q[2]^2-q[3]^2+q[4]^2]
end

function r2q(r::Matrix{<:Number})::Vector{Float64}
    c = 0.5*(tr(r)-1)
    if c == 1
        return [0, 0, 0, 1]
    elseif c == -1
        m = r + I
        for i in 1:3
            v = m[:, i]
            if norm(v) > eps()
                return 1/norm(v)*vcat([0], v)
            end
        end
    end
    q = 1/(2*sqrt(tr(r) + 1))*[tr(r)+1, r[3, 2]-r[2, 3], r[1, 3] - r[3, 1], r[2, 1] - r[1, 2]]
    return 1/norm(q)*q
end


function rat_approx(n::Float64, tol::Float64)::Rational{BigInt}
    function fexp(f)
        return f != 0 ? floor(log10(abs(f))) : 0
    end
    return tol >= 1 ? Int(floor(n)) : Int(floor(n*10^(-fexp(tol))))//Int(10^(-fexp(tol)))
end

n = randn()
tol = 1e-5
rat = rat_approx(n, tol)
abs(rat - n) < tol




function exact_cs(angle::Float64, tol::Float64)::Tuple{Rational{BigInt}, Rational{BigInt}}
    t = rat_approx(tan(angle/2), tol/10)
    return (1-t^2)/(1+t^2), 2*t/(1+t^2)
end

angle = randn()
tol = 1e-5
c, s = exact_cs(angle, tol)
c^2 + s^2 == 1
abs(atan(s, c) - angle) < tol




function exact_rot(q::Vector{Float64}, tol::Float64)::Matrix{Rational{BigInt}}
    tol_q = tol
    while true
        q_rat = [rat_approx(n, tol_q) for n in q]
        r = q2r(q_rat)
        if norm(r - q2r(q)) < tol
            return r
        else
            tol_q /= 10
        end
    end
end

q = randn(4)
tol = 1e-5
r = exact_rot(q, tol)
r'*r == I
det(r) == 1



function rational_mechanism(mechanism::Dict{String, Vector{Float64}}, tol::Float64)::Dict{String, Vector{Real}}
    rat_mechanism = Dict{String, Vector{Real}}("theta offset" => zeros(6),
                                                  "d" => zeros(6),
                                                  "a" => zeros(6),
                                                  "cos alpha" => zeros(6),
                                                  "sin alpha" => zeros(6))
    for i in 1:6
        rat_mechanism["theta offset"][i] = mechanism["theta offset"][i]
        rat_mechanism["d"][i] = rat_approx(mechanism["d"][i], tol)
        rat_mechanism["a"][i] = rat_approx(mechanism["a"][i], tol)
        rat_mechanism["cos alpha"][i], rat_mechanism["sin alpha"][i] = exact_cs(mechanism["alpha"][i], tol)
    end
    return rat_mechanism
end

function rational_pose(pose::Dict{String, Vector{Float64}}, tol::Float64)::Matrix{Rational{BigInt}}
    t = [rat_approx(n, tol) for n in pose["t"]]
    r = exact_rot(pose["q"], tol)
    e = [r t; 0 0 0 1]
end

pose = Dict{String, Vector{Float64}}("q" => randn(4), "t" => randn(3))
tol = 1e-5
rat_pose = rational_pose(pose, tol)
r = rat_pose[1:3,1:3]
t = rat_pose[1:3,4]
r'*r == I
det(r) == 1
norm(r - q2r(pose["q"])) < tol
norm(t - pose["t"]) < 10*tol


function rat2nQ!(M::Matrix)
    for i in 1:size(M, 1)
        for j in 1:size(M, 2)
            if typeof(M[i, j]) != spoly{n_Q}
                M[i, j] = QQ(M[i, j])
            end
        end
    end
end

function dh_matrix_sym(cos_th, sin_th, d::Number, a::Number, cos_alpha::Number, sin_alpha::Number)::Matrix
    t1 = [cos_th -sin_th 0 0; sin_th cos_th 0 0; 0 0 1 d; 0 0 0 1]
    t2 = Matrix{Any}(sparsify([1 0 0 a; 0 cos_alpha -sin_alpha 0; 0 sin_alpha cos_alpha 0; 0 0 0 1], 1e-5))
    # if typeof(cos_th) == spoly{n_Q}
    #     rat2nQ!(t1)
    #     rat2nQ!(t2)
    # end
    return t1 * t2
end

function inv_dh_matrix_sym(cos_th, sin_th, d::Number, a::Number, cos_alpha::Number, sin_alpha::Number)::Matrix
    t1 = [cos_th sin_th 0 0; -sin_th cos_th 0 0; 0 0 1 -d; 0 0 0 1]
    t2 = Matrix{Any}(sparsify([1 0 0 -a; 0 cos_alpha sin_alpha 0; 0 -sin_alpha cos_alpha 0; 0 0 0 1], 1e-5))
    if typeof(cos_th) == spoly{n_Q}
        rat2nQ!(t1)
        rat2nQ!(t2)
    end
    return t2 * t1
end

function transformation_sym(mechanism::Dict{String, <:Vector{<:Number}}, idx_from::Int64, idx_to::Int64, c::Vector, s::Vector)::Matrix
    # symb = typeof(c[1]) == spoly{n_Q} ? true : false
    symb = false
    t = I
    for i in idx_to+1:idx_from
        cos_th, sin_th = c[i], s[i]
        d = mechanism["d"][i]
        a = mechanism["a"][i]
        if symb
            cos_alpha, sin_alpha = mechanism["cos alpha"][i], mechanism["sin alpha"][i]
            t = t * dh_matrix_sym(cos_th, sin_th, d, a, cos_alpha, sin_alpha)
        else
            alpha = mechanism["alpha"][i]
            t = t * dh_matrix_sym(cos_th, sin_th, d, a, cos(alpha), sin(alpha))
        end
    end
    return t
end

function inv_transformation_sym(mechanism::Dict{String, <:Vector{<:Number}}, idx_from::Int64, idx_to::Int64, c::Vector, s::Vector)::Matrix
    rat = typeof(mechanism["d"][1]) <: Rational ? true : false
    t = I
    for i in idx_from:-1:idx_to+1
        cos_th, sin_th = c[i], s[i]
        d = mechanism["d"][i]
        a = mechanism["a"][i]
        if rat
            cos_alpha, sin_alpha = mechanism["cos alpha"][i], mechanism["sin alpha"][i]
            t = t * inv_dh_matrix_sym(cos_th, sin_th, d, a, cos_alpha, sin_alpha)
        else
            alpha = mechanism["alpha"][i]
            t = t * inv_dh_matrix_sym(cos_th, sin_th, d, a, cos(alpha), sin(alpha))
        end
    end
    return t
end

function dh_matrix(theta::Float64, d::Float64, a::Float64, alpha::Float64)::Matrix{Float64}
    t1 = sparsify([cos(theta) -sin(theta) 0 0; sin(theta) cos(theta) 0 0; 0 0 1 d; 0 0 0 1], 1e-5)
    t2 = sparsify([1 0 0 a; 0 cos(alpha) -sin(alpha) 0; 0 sin(alpha) cos(alpha) 0; 0 0 0 1], 1e-5)
    return t1 * t2
end

function fkt(mechanism::Dict{String, Vector{Float64}}, joints::Vector{Float64})::Dict{String, Vector{Float64}}
    t = I
    for i in 1:6
        theta = mechanism["theta offset"][i] + joints[i]
        d = mechanism["d"][i]
        a = mechanism["a"][i]
        alpha = mechanism["alpha"][i]
        t = t * dh_matrix(theta, d, a, alpha)
    end
    return Dict{}("q" => r2q(t[1:3, 1:3]), "t" => t[1:3, 4])
end

function random_mechanism()
    return Dict{String, Vector{Float64}}("theta offset" => randn(6),
                                         "d" => randn(6),
                                         "a" => randn(6),
                                         "alpha" => randn(6))
end
