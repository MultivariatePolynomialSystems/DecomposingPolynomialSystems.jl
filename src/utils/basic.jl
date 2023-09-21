using HomotopyContinuation: Expression, Variable
export a2p, p2a, M2V, M2VV, M2VM, V2M, V2Mt, VV2M, VM2M, VM2V, xx, xx2v, eye

function a2p(M)
    return [M; ones(eltype(M), 1, size(M, 2))]
end

function p2a(M)
    Ma = M./reshape(M[end,:], 1, size(M, 2))
    return Ma[1:end-1,:]
end

function M2V(M::Array)
    return M[:]
end

function M2VV(M)
    return [M[:,i] for i in axes(M, 2)]
end

function M2VM(M::Matrix)
    return [reshape(M[:,i], size(M,1), 1) for i in axes(M, 2)]
end

function V2M(V::Vector)
    return reshape(V, :, 1)
end

function V2Mt(V::Vector)
    return hcat(V...)
end

function VV2M(V::Vector)
    return hcat(V...)
end

function VM2M(V::Vector)
    return hcat(V...)
end

function VM2V(V::Vector)
    return M2V(hcat(V...))
end

function xx(v)
    return [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]
end

function xx2v(xx)
    return [-xx[2,3], xx[1,3], -xx[1,2]]
end

function eye(T, n)
    return Matrix{T}(I(n))
end