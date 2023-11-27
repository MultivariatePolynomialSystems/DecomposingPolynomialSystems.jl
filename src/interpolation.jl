export rational_function,
    polynomial_function,
    remove_zero_nums_and_denoms

using LinearAlgebra: norm, dot

function rational_function(
    coeffs::AbstractVector{<:Number},
    num_mons::MonomialVector,
    denom_mons::MonomialVector;
    logging::Bool=false,
    tol::Float64=1e-5
)
    n_num_mons, n_denom_mons = length(num_mons), length(denom_mons)
    @assert length(coeffs) == n_num_mons + n_denom_mons
    num, denom = coeffs[1:n_num_mons], coeffs[n_num_mons+1:end]
    if norm(num) < tol
        @warn "Numerator close to zero"
    end
    @assert norm(denom) > tol
    nonzero_ids = findall(!iszero, denom)
    nonzero_denom = denom[nonzero_ids]
    if all(x->x==nonzero_denom[1], nonzero_denom)
        num /= nonzero_denom[1]
        sparsify!(num, tol; digits=1)
        denom[nonzero_ids] = ones(length(nonzero_ids))
    end
    num, denom = simplify_numbers(num), simplify_numbers(denom)
    if logging
        println("numerator = ", sum(to_expressions(num_mons).*num))
        println("denominator = ", sum(to_expressions(denom_mons).*denom))
    end
    p = sum(to_expressions(num_mons).*num)
    q = sum(to_expressions(denom_mons).*denom)
    logging && println("rational function = ", p/q)
    return p/q
end

function polynomial_function(
    coeffs::AbstractVector{<:Number},
    mons::MonomialVector;
    logging::Bool=false
)
    @assert length(coeffs) == length(mons.mds)
    p = sum(to_expressions(mons).*coeffs)
    logging && println("polynomial = ", p)
    return p
end

function check_func_type(func_type::String)
    if (func_type != "polynomial" && func_type != "rational")
        error("func_type argument must be either \"polynomial\" or \"rational\"")
    end
end

function reconstruct_function(
    coeffs::AbstractVector{<:Number},
    mons::MonomialVector;
    func_type::String,
    logging::Bool=true,
    tol::Float64=1e-5
)
    check_func_type(func_type)
    if func_type == "rational"
        return rational_function(coeffs, mons, mons, logging=logging, tol=tol)
    else
        return polynomial_function(coeffs, mons, logging=logging)
    end
end

function remove_zero_nums_and_denoms(
    coeffs::AbstractMatrix{<:Number},
    num_mons::MonomialVector,
    denom_mons::MonomialVector;
    logging::Bool=false
)

    reasonable_rows = []
    n_num_mons, n_denom_mons = length(num_mons.mds), length(denom_mons.mds)
    @assert size(coeffs, 2) == n_num_mons + n_denom_mons
    for i in axes(coeffs, 1)
        if (!iszero(coeffs[i, 1:n_num_mons]) && !iszero(coeffs[i, n_num_mons+1:end]))
            push!(reasonable_rows, i)
        elseif logging
            if iszero(coeffs[i, 1:n_num_mons])
                println(
                    "Denominator removed: ",
                    polynomial_function(coeffs[i, n_num_mons+1:end], denom_mons)
                )
            else
                println(
                    "Numerator removed: ",
                    polynomial_function(coeffs[i, 1:n_num_mons], num_mons)
                )
            end
        end
    end
    return coeffs[reasonable_rows, :]
end

function remove_zero_nums_and_denoms(
    coeffs::AbstractMatrix{<:Number},
    mons::MonomialVector;
    logging::Bool=false
)
    return _remove_zero_nums_and_denoms(coeffs, mons, mons, logging=logging)
end

function rational_functions(
    coeffs::AbstractMatrix{<:Number},
    num_mons::MonomialVector,
    denom_mons::MonomialVector;
    logging::Bool=true,
    tol::Float64=1e-5
)
    for k in eachindex(axes(coeffs, 1))
        rational_function(coeffs[k,:], num_mons, denom_mons, logging=logging, tol=tol)
        println()
    end
end

function reconstruct_functions(
    coeffs::AbstractMatrix{<:Number},
    mons::MonomialVector;
    func_type::String,
    logging::Bool=true,
    tol::Float64=1e-5
)
    
end

function good_representative(coeffs::AbstractMatrix{<:Number})
    return coeffs[findmin(vec(sum(coeffs .!= 0, dims=2)))[2], :]
end

# NOT READY YET...
function best_representative(rational_functions::AbstractMatrix{<:Number}, tol::Float64)
    n_mons = Int(size(rational_functions, 2)/2)
    A = rational_functions[:,n_mons+1:end]
    A = [-1 zeros(1, n_mons-1); A]
    R = Matrix{CC}(transpose(nullspace(transpose(A))))
    R = sparsify(rref(R, tol), tol, digits=1) # 1*R[1,:] + r2*R[2,:] + ... are the solutions
    a = R[:,2:end]*rational_functions[:,1:n_mons] # 1*a[1,:] + r2*a[2,:] + ... are possible numerators
    return nothing
end

function vandermonde_matrix(
    values::AbstractVector{<:Number},
    eval_mons::AbstractMatrix{<:Number},
    func_type::String
)
    check_func_type(func_type)
    if func_type == "polynomial"
        return [transpose(eval_mons) -values]
    elseif func_type == "rational"
        return [transpose(eval_mons) -transpose(eval_mons).*values]
    end
end

# supposes each md in mds is a multidegree in both unknowns and parameters
# TODO: The implementation below (with _) is more efficient (approx 2x),
# TODO: since it exploits the sparsity of multidegrees. REMOVE THIS METHOD?
# function HomotopyContinuation.evaluate(mons::MonomialVector, samples::Samples)
#     solutions = samples.solutions
#     parameters = samples.parameters

#     n_unknowns, n_sols, n_instances = size(solutions)
#     mds = mons.mds
#     n_mds = length(mds)

#     evaluated_mons = zeros(CC, n_mds, n_sols, n_instances)
#     for i in 1:n_instances
#         params = view(parameters, :, i)
#         params_eval = [prod(params.^md[n_unknowns+1:end]) for md in mds]
#         sols = view(solutions, :, :, i)
#         for j in 1:n_mds
#             evaluated_mons[j, :, i] = vec(prod(sols.^mds[j][1:n_unknowns], dims=1)).*params_eval[j]
#         end
#     end
#     return evaluated_mons
# end

function interpolate(
    A::AbstractMatrix{CC},
    mons::MonomialVector;
    func_type::String,
    tol::Float64=1e-5
)
    check_func_type(func_type)
    @assert size(A, 1) >= size(A, 2)
    if func_type == "polynomial"
        @assert size(A, 2) == length(mons)
    else
        @assert size(A, 2) == 2*length(mons)
    end

    C = Matrix{CC}(transpose(nullspace(A)))
    C = sparsify(rref(C, tol), tol, digits=1)
    if size(C, 1) == 0
        return nothing
    elseif size(C, 1) == 1
        return reconstruct_function(C[1, :], mons, func_type=func_type)
    else
        # choose good/best representative
        return reconstruct_function(C[1, :], mons, func_type=func_type)
    end
end

# If samples are generated from a proper subvariety X ⊂ Cⁿ, then the method returns just one representative f
# of the equivalence class of functions f+I(X) with the given values. Obviously, if X = Cⁿ, then there is only
# one such function.
function interpolate(
    values::AbstractVector{CC},
    mons::MonomialVector,
    eval_mons::AbstractMatrix{CC};
    func_type::String,
    tol::Float64=1e-5
)
    check_func_type(func_type)
    A = vandermonde_matrix(values, eval_mons, func_type)
    return interpolate(A, mons, func_type=func_type, tol=tol)
end

function interpolate_vanishing_polynomials(
    samples::Samples,
    vars::Vector{Variable},
    mons::MonomialVector;
    tol::Float64=1e-5
)
    A = evaluate_monomials_at_samples(mons, samples, vars)
    A = reshape(A, size(A, 1), size(A, 2)*size(A, 3))
    N = Matrix{CC}(transpose(nullspace(transpose(A))))
    N = sparsify(rref(N, tol), tol, digits=1)
    return N*mons
end

function interpolate_vanishing_polynomials(
    samples::AbstractMatrix{CC},
    vars::Vector{Variable},
    mons::MonomialVector;
    tol::Float64=1e-5
)
    @assert size(samples, 1) == length(vars)
    @assert size(samples, 2) >= length(mons)
    A = evaluate_monomials_at_samples(mons, samples, vars)
    N = Matrix{CC}(transpose(nullspace(transpose(A))))
    N = sparsify(rref(N, tol), tol, digits=1)
    return N*mons
end

function interpolate_dense(
    samples::AbstractMatrix{CC},
    values::AbstractVector{CC};
    vars::Vector{Variable},
    degree::Int,
    func_type::String,
    tol::Float64=1e-5
)
    
end

function interpolate_sparse(
    samples::AbstractMatrix{CC},
    values::AbstractVector{CC};
    vars::Vector{Variable},
    degree::Int,
    func_type::String,
    tol::Float64=1e-5
)
    
end

function interpolate(
    samples::AbstractMatrix{CC},
    values::AbstractVector{CC};
    vars::Vector{Variable},
    degree::Int,
    func_type::String,
    method::String,
    tol::Float64=1e-5
)
    if method == "dense"
        return interpolate_dense(samples, values; vars=vars, degree=degree, func_type=func_type, tol=tol)
    elseif method == "sparse"
        return interpolate_sparse(samples, values; vars=vars, degree=degree, func_type=func_type, tol=tol)
    else
        error("Not supported method '", method, "' for interpolation")
    end
end







#----------- SPARSE INTERPOLATION TESTING -------------

function inv_chinese(d::Int, p::Vector{Int})
    n = length(p)
    m = prod(p)
    D = zeros(Int, n)
    for i in 1:n
        D[i] = mod(d*invmod(Int(m/p[i]), p[i]), p[i])
    end
    return D
end

function companion_matrix(λ::Vector{CC})::Matrix{CC}
    d = length(λ)
    return [[zeros(1,d-1); eye(d-1)] -λ]
end

function f(x, y)
    return x^20*y + 2*x^2 + 3*x + 4*y^11 + 5
end

function interpolate_sparse()

end

# t = 5 # number of terms in f
# p₁, p₂ = 139, 193 # pᵢ > deg(fₓᵢ)
# m = p₁*p₂
# ω₁, ω₂ = cispi(2/p₁), cispi(2/p₂) # pᵢ-th roots of unity
# ω = cispi(2/m)

# α = [f(ω₁^s, ω₂^s) for s in 0:2*t-1]

# H₀ = VV2M([[α[j] for j in i:i+t-1] for i in 1:t])
# h = -[α[i] for i in t+1:2*t]
# Hh = [H₀ -h]
# N = nullspace(Hh)
# λ = M2V(p2a(N))

# M = companion_matrix(λ) 
# b = eigvals(M)

# d = [mod(Int(round(log(ω, b[i]))), m) for i in 1:length(b)]
# e = [inv_chinese(d[i], [p₁, p₂]) for i in 1:length(d)]

# B = [transpose(VV2M([b.^i for i in 0:t-1])) -α[1:t]]
# N = nullspace(B)
# c = round.(M2V(p2a(N)))

# @var x y
# f_approx = dot(c, [prod([x,y].^e[i]) for i in 1:length(e)])
