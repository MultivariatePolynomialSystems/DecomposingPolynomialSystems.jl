export rref

function rref!(A::Matrix{T}, ɛ=T <: Union{Rational,Integer} ? 0 : eps(norm(A,Inf))) where T
    nr, nc = size(A)
    i = j = 1
    while i <= nr && j <= nc
        (m, mi) = findmax(abs.(A[i:nr,j]))
        mi = mi+i - 1
        if m <= ɛ
            if ɛ > 0
                A[i:nr,j] .= zero(T)
            end
            j += 1
        else
            for k=j:nc # swap i-th and mi-th rows
                A[i, k], A[mi, k] = A[mi, k], A[i, k]
            end
            d = A[i,j] # == m
            for k = j:nc
                A[i,k] /= d # make pivot = 1
            end
            for k = 1:nr # subtract multiple of pivot row from every row except i-th
                if k != i
                    d = A[k,j]
                    for l = j:nc
                        A[k,l] -= d*A[i,l]
                    end
                end
            end
            i += 1
            j += 1
        end
    end
    A
end

rrefconv(::Type{T}, A::Matrix) where {T} = rref!(copyto!(similar(A, T), A))

rref(A::Matrix{T}) where {T} = rref!(copy(A))
rref(A::Matrix{T}) where {T <: Complex} = rrefconv(CC, A)
rref(A::Matrix{CC}, tol::Float64) = rref!(copy(A), tol)
rref(A::Matrix{<:Number}, tol::Float64) = rref!(copy(A), tol)
rref(A::Matrix{T}) where {T <: Union{Integer, Float16, Float32}} = rrefconv(Float64, A)
rref(A::AbstractMatrix) = rref(Matrix(A))
