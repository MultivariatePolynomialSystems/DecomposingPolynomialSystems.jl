export sparsify!,
    simplify_numbers,
    eye, a2p, p2a,
    num_mons, num_mons_upto

a2p(M::AbstractMatrix{<:Number}) = [M; ones(eltype(M), 1, size(M, 2))]
p2a(M::AbstractMatrix{<:Number}) = (M./M[end:end,:])[1:end-1,:]

M2VV(M::AbstractMatrix) = [M[:,i] for i in axes(M, 2)]
M2VM(M::AbstractMatrix) = [reshape(M[:,i], size(M,1), 1) for i in axes(M, 2)]

function Base.copyto!(M::AbstractMatrix{T}, v::AbstractVector{AbstractVector{T}}; dim::Integer) where {T}
    for i in eachindex(v)
        copyto!(selectdim(M, dim, i), v[i])
    end
end

xx(v) = [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]
xx2v(xx) = [-xx[2,3], xx[1,3], -xx[1,2]]
eye(T, n::Integer) = Matrix{T}(I(n))

prodpow(v::AbstractVector, e::AbstractSparseVector) = prod(v[e.nzind].^e.nzval)

num_mons(n::Integer, d::Integer) = n > 0 ? binomial(Int(n - 1 + d), Int(d)) : 0
num_mons_upto(n::Integer, d::Integer) = n > 0 ? binomial(Int(n + d), Int(d)) : 0

# TODO: test this
function sparsify!(v::AbstractVector{<:Number}, tol::Real; digits::Integer=0)
    for j in eachindex(v)
        if abs(imag(v[j])) < tol
            v[j] = real(v[j])
        elseif abs(round(imag(v[j]); digits=digits) - imag(v[j])) < tol
            v[j] = real(v[j]) + round(imag(v[j]); digits=digits)*im
        end
        if abs(real(v[j])) < tol
            v[j] = imag(v[j])*im
        elseif abs(round(real(v[j]); digits=digits) - real(v[j])) < tol
            v[j] = round(real(v[j]); digits=digits) + imag(v[j])*im
        end
    end
    return v
end

function sparsify!(M::AbstractMatrix{<:Number}, tol::Real; digits::Integer=0)
    for r in eachrow(M)
        sparsify!(r, tol; digits=digits)
    end
    return M
end

function simplify_numbers(v::Vector{<:Number})
    v = Vector{Number}(v)
    for (i, vᵢ) in enumerate(v)
        try
            v[i] = Integer(vᵢ)
        catch
            try
                v[i] = Real(vᵢ)
            catch
                try
                    v[i] = Complex{Integer}(vᵢ)
                catch
                end
            end
        end
    end
    return v
end

function to_ordinal(n::Integer)::String
    if mod(n, 10) == 1
        mod(n, 100) == 11 && return "$(n)th"
        return "$(n)st"
    end
    if mod(n, 10) == 2
        mod(n, 100) == 12 && return "$(n)th"
        return "$(n)nd"
    end
    if mod(n, 10) == 3
        mod(n, 100) == 13 && return "$(n)th"
        return "$(n)rd"
    end
    return "$(n)th"
end

function subscript(n::Integer)::String
    c = n < 0 ? [Char(0x208B)] : []
    for d in reverse(digits(abs(n)))
        push!(c, Char(0x2080+d))
    end
    return join(c)
end

function superscript(n::Integer)::String
    c = n < 0 ? [Char(0x207B)] : []
    for d in reverse(digits(abs(n)))
        if d == 0 push!(c, Char(0x2070)) end
        if d == 1 push!(c, Char(0x00B9)) end
        if d == 2 push!(c, Char(0x00B2)) end
        if d == 3 push!(c, Char(0x00B3)) end
        if d > 3 push!(c, Char(0x2070+d)) end
    end
    return join(c)
end

# TODO: eachcol?
Base.findfirst(
    v::AbstractVector{<:Number},
    M::AbstractMatrix{<:Number};
    tol::Real=1e-5
) = findfirst(i->norm(M[:,i]-v)<tol, axes(M,2))

take_rows(
    f::Function,
    M::AbstractMatrix{T}
) where {T} = M[[f(r) for r in eachrow(M)], :]

function column_diffs(M::AbstractMatrix{T}) where {T<:Number}
    M = M - M[:,1]*ones(T, 1, size(M,2))
    return M[:,2:end]
end

phrase(i::Integer, word::String) = i == 1 ? "$(i) $(word)" : "$(i) $(word)s"

# function exprDet(M; expnd=true)
#     n = size(M, 1)
#     @assert(n == size(M, 2))
#     if n == 1 return M[1,1]
#     elseif n == 2
#         if expnd return expand(M[1,1]*M[2,2]-M[1,2]*M[2,1])
#         else return M[1,1]*M[2,2]-M[1,2]*M[2,1]
#         end
#     else
#         if expnd return expand(sum([(-1)^(j+1)*M[1,j]*exprDet(M[2:n,filter(k->k!=j,convert(Vector,1:n))], expnd=expnd) for j=1:n]))
#         else return sum([(-1)^(j+1)*M[j,1]*exprDet(M[filter(k->k!=j,convert(Vector,1:n)), 2:n], expnd=expnd) for j=1:n])
#         end
#     end
# end

# function exprAdj(M::Matrix{Expression}; expnd=true)
#     m, n = size(M)
#     @assert m == n
#     A = Matrix{Expression}(zeros(n, n))
#     for i=1:n for j=1:n A[i,j] = (-1)^(i+j)*exprDet(M[1:end .!= i, 1:end .!= j], expnd=expnd) end end
#     return Matrix{Expression}(transpose(A))
# end

# function exprInv(M::Matrix{Expression})
#     m, n = size(M)
#     @assert m == n
#     return 1/exprDet(M)*exprAdj(M)
# end


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


# # Returns list of monomials up to degree d
# function get_monomials(x::Vector{Variable}, d::Int)::Vector{Expression}
#     mons = Vector{Expression}([1])
#     n_vars = length(x)
#     append!(mons, x)
#     k = 1
#     for i = 1:d-1
#         M = mons[(k + 1):(k + num_mons(n_vars, i))] * transpose(x)
#         for j = 1:n_vars
#             append!(mons, M[1:num_mons(j, i), j])
#         end
#         k += num_mons(n_vars, i)
#     end
#     return mons
# end

# function next_deg_mons(x::Vector{Variable}, mons::Vector{Expression}, d::Int)::Vector{Expression}
#     new_mons = Vector{Expression}([])
#     n_vars = length(x)
#     M = mons[(end-num_mons(n_vars, d-1)+1):end] * transpose(x)
#     for j = 1:n_vars
#         append!(new_mons, M[1:num_mons(j, d-1), j])
#     end
#     return new_mons
# end

# function get_monomials_fixed_degree(x::Vector{Variable}, d::Int)::Vector{Expression}
#     if d == 0 return [1] end
#     if d == 1 return x end
#     mons = Vector{Expression}([1])
#     n_vars = length(x)
#     append!(mons, x)
#     k = 1
#     for i = 1:d-1
#         M = mons[(k + 1):(k + num_mons(n_vars, i))] * transpose(x)
#         for j = 1:n_vars
#             append!(mons, M[1:num_mons(j, i), j])
#         end
#         k += num_mons(n_vars, i)
#         if i == d-1
#             return mons[(k+1):(k + num_mons(n_vars, d))]
#         end
#     end
# end

# # Returns list of monomials up to degree d for factorizing map
# function get_monomials_factorization(x::Vector{Variable}, d::Int, n_params::Int)::Vector{Expression}
#     mons = Vector{Expression}([1])
#     n_vars = length(x)
#     n_unknowns = n_vars - n_params
#     append!(mons, x[1:n_unknowns])
#     k = 1
#     for i = 1:d-1
#         n_old_mons = num_mons(n_vars, i) - num_mons(n_params, i)
#         M = mons[(k + 1):(k + n_old_mons)] * transpose(x)
#         for j = 1:n_vars
#             n_new_mons = min(num_mons(j, i) - num_mons(j-n_unknowns, i), size(M, 1))
#             append!(mons, M[1:n_new_mons, j])
#         end
#         k += n_old_mons
#     end
#     popfirst!(mons)
#     return mons
# end


# # rewrite
# # TODO: what is this?
# function v2SymM(v)
#     n = Int((-1+sqrt(1+8*length(v)))/2)
#     M = zeros(n, n)
#     k = 1
#     for i=1:n
#         for j=i:n
#             M[i,j] = M[j,i] = v[k]
#             k += 1
#         end
#     end
#     return M
# end