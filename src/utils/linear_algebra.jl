using HomotopyContinuation: Expression, Variable
export q2R, c2R, exprDet, sparsify, num_mons, num_mons_upto, get_monomials, get_monomials_fixed_degree

function exprDet(M; expnd=true)
    n = size(M, 1)
    @assert(n == size(M, 2))
    if n == 1 return M[1,1]
    elseif n == 2
        if expnd return expand(M[1,1]*M[2,2]-M[1,2]*M[2,1])
        else return M[1,1]*M[2,2]-M[1,2]*M[2,1]
        end
    else
        if expnd return expand(sum([(-1)^(j+1)*M[1,j]*exprDet(M[2:n,filter(k->k!=j,convert(Vector,1:n))], expnd=expnd) for j=1:n]))
        else return sum([(-1)^(j+1)*M[j,1]*exprDet(M[filter(k->k!=j,convert(Vector,1:n)), 2:n], expnd=expnd) for j=1:n])
        end
    end
end

function exprAdj(M::Matrix{Expression}; expnd=true)
    m, n = size(M)
    @assert m == n
    A = Matrix{Expression}(zeros(n, n))
    for i=1:n for j=1:n A[i,j] = (-1)^(i+j)*exprDet(M[1:end .!= i, 1:end .!= j], expnd=expnd) end end
    return Matrix{Expression}(transpose(A))
end

function exprInv(M::Matrix{Expression})
    m, n = size(M)
    @assert m == n
    return 1/exprDet(M)*exprAdj(M)
end

function q2R(q)
    a,b,c,d = q
    return 1/(a^2+b^2+c^2+d^2)*[a^2+b^2-c^2-d^2 2*(b*c-a*d) 2*(b*d+a*c); 2*(b*c+a*d) a^2-b^2+c^2-d^2 2*(c*d-a*b); 2*(b*d-a*c) 2*(c*d+a*b) a^2-b^2-c^2+d^2]
end

function q2c(q)
    return -[q[2]/q[1], q[3]/q[1], q[4]/q[1]]
end

function c2R(c)
    x, y, z = c
    return 1/(1+x^2+y^2+z^2)*[1+x^2-y^2-z^2 2*(x*y-z) 2*(y+x*z); 2*(x*y+z) 1-x^2+y^2-z^2 2*(y*z-x); 2*(x*z-y) 2*(x+y*z) 1-x^2-y^2+z^2]
end

function R2c(R)
    return 1/(1+tr(R))*xx2v((R'-R))
end

function ct2P(c,t)
    return [c2R(c) t]
end

function ct2Prad(c,t)
    R = c2R(c)
    return [R[1:2,:] t]
end

function c2sR(c)
    x, y, z = c
    return [1+x^2-y^2-z^2 2*(x*y-z) 2*(y+x*z); 2*(x*y+z) 1-x^2+y^2-z^2 2*(y*z-x); 2*(x*z-y) 2*(x+y*z) 1-x^2-y^2+z^2]
end

# function c2R(c)
#     return (I - xx(c))*inv(I + xx(c))
# end

function ct2sP(c,t)
    sR = c2sR(c)
    return [sR (1+c[1]^2+c[2]^2+c[3]^2)*t]
end

function ct2sPrad(c,t)
    sR = c2sR(c)
    return [sR[1:2,:] (1+c[1]^2+c[2]^2+c[3]^2)*t]
end

# rewrite
function v2SymM(v)
    n = Int((-1+sqrt(1+8*length(v)))/2)
    M = zeros(n, n)
    k = 1
    for i=1:n
        for j=i:n
            M[i,j] = M[j,i] = v[k]
            k += 1
        end
    end
    return M
end

function sparsify(x::Vector{<:Number}, tol::Float64; digits::Int64=0)
    for j in 1:length(x)
        if abs(imag(x[j])) < tol
            x[j] = real(x[j])
        elseif abs(round(imag(x[j]), digits=digits) - imag(x[j])) < tol
            x[j] = real(x[j]) + round(imag(x[j]), digits=digits)*im
        end
        if abs(real(x[j])) < tol
            x[j] = imag(x[j])*im
        elseif abs(round(real(x[j]), digits=digits) - real(x[j])) < tol
            x[j] = round(real(x[j]), digits=digits) + imag(x[j])*im
        end
    end
    return x
end

function sparsify(A::Matrix{<:Number}, tol::Float64; digits::Int64=0)
    for i in 1:size(A, 2)
        x = sparsify(A[:, i], tol, digits=digits)
        A[:, i] = x
    end
    return A
end

# number of monomials of degree d in n variables
function num_mons(n, d)
    return n > 0 ? binomial(n - 1 + d, d) : 0
end

function num_mons_upto(n, d)
    return n > 0 ? binomial(n + d, d) : 0
end

# Returns list of monomials up to degree d
function get_monomials(x::Vector{Variable}, d::Int)::Vector{Expression}
    mons = Vector{Expression}([1])
    n_vars = length(x)
    append!(mons, x)
    k = 1
    for i = 1:d-1
        M = mons[(k + 1):(k + num_mons(n_vars, i))] * transpose(x)
        for j = 1:n_vars
            append!(mons, M[1:num_mons(j, i), j])
        end
        k += num_mons(n_vars, i)
    end
    return mons
end

function next_deg_mons(x::Vector{Variable}, mons::Vector{Expression}, d::Int)::Vector{Expression}
    new_mons = Vector{Expression}([])
    n_vars = length(x)
    M = mons[(end-num_mons(n_vars, d-1)+1):end] * transpose(x)
    for j = 1:n_vars
        append!(new_mons, M[1:num_mons(j, d-1), j])
    end
    return new_mons
end

function get_monomials_fixed_degree(x::Vector{Variable}, d::Int)::Vector{Expression}
    if d == 0 return [1] end
    if d == 1 return x end
    mons = Vector{Expression}([1])
    n_vars = length(x)
    append!(mons, x)
    k = 1
    for i = 1:d-1
        M = mons[(k + 1):(k + num_mons(n_vars, i))] * transpose(x)
        for j = 1:n_vars
            append!(mons, M[1:num_mons(j, i), j])
        end
        k += num_mons(n_vars, i)
        if i == d-1
            return mons[(k+1):(k + num_mons(n_vars, d))]
        end
    end
end

# Returns list of monomials up to degree d for factorizing map
function get_monomials_factorization(x::Vector{Variable}, d::Int, n_params::Int)::Vector{Expression}
    mons = Vector{Expression}([1])
    n_vars = length(x)
    n_unknowns = n_vars - n_params
    append!(mons, x[1:n_unknowns])
    k = 1
    for i = 1:d-1
        n_old_mons = num_mons(n_vars, i) - num_mons(n_params, i)
        M = mons[(k + 1):(k + n_old_mons)] * transpose(x)
        for j = 1:n_vars
            n_new_mons = min(num_mons(j, i) - num_mons(j-n_unknowns, i), size(M, 1))
            append!(mons, M[1:n_new_mons, j])
        end
        k += n_old_mons
    end
    popfirst!(mons)
    return mons
end

# function get_specific_monomials(x::Vector{Variable}, d::Int64)::Vector{Expression}
#     mons = Vector{Expression}([1])
#     n_vars = length(x)
#     append!(mons, x)
#     n_appended = n_vars
#     for i = 1:d-1
#         M = mons[end-n_appended+1:end] * transpose(x)
#         for j = i+1:n_vars
#             append!(mons, M[1:j-i, j])
#         end
#         n_appended = Int((n_vars-i)*(n_vars-i+1)/2)
#     end
#     return mons
# end

# function get_monomials(x::Vector{Variable}, d::Int64)::Vector{Expression}
#     return [prod(x.^exp) for exp in collect(multiexponents(length(x), d))]
# end
#
# function get_monomials(x::Vector{Variable}, degrees::Vector{Int64})::Vector{Expression}
#     mons = []
#     n = length(x)
#     for d in degrees
#         append!(mons, get_monomials(x, d))
#     end
#     return mons
# end
