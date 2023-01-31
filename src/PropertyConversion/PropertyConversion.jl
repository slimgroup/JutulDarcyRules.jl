export TransToK, KtoTrans, K1to3

K1to3(Kx::AbstractArray{T}) where T = vcat(vec(Kx)', vec(Kx)', vec(Kx)')

function TransToK(g::CartesianMesh, Trans::AbstractVector{T}) where T<:Number

    # ## Set up L-BFGS
    function fg(F, G, x)
        grads = gradient(Flux.params(x)) do
            F = norm(KtoTrans(g, x)-Trans)^2 / norm(Trans)^2
        end
        G .= grads[x]
        return F
    end

    Kx = vcat([vec(g.deltas[i]^2 / prod(g.deltas) * mean(Trans) * ones(prod(g.dims)))' for i = 1:3]...)
    result = optimize(Optim.only_fg!(fg), Kx, LBFGS(), Optim.Options(f_tol=T(0)))
    Kxinv = result.minimizer
    return Kxinv

end

function KtoTrans(g::CartesianMesh, K::AbstractArray{T}) where {T<:Number}
    d = g.deltas
    n = g.dims
    return vcat(vec(KtoTransx(reshape(K[1,:], n), d)),
                vec(KtoTransy(reshape(K[2,:], n), d)),
                vec(KtoTransz(reshape(K[3,:], n), d)))
end

function KtoTransx(K::AbstractArray{T, 3}, d::NTuple{3, T}) where {T<:Number}
    Kreci = T(1) ./ K
    return T(2) * prod(d) / d[1]^2 ./ (Kreci[1:end-1,:,:] + Kreci[2:end,:,:])
end

function KtoTransy(K::AbstractArray{T, 3}, d::NTuple{3, T}) where {T<:Number}
    Kreci = T(1) ./ K
    return permutedims(T(2) .* prod(d) / d[2]^2 ./ (Kreci[:,1:end-1,:] + Kreci[:,2:end,:]), [1,3,2])
end

function KtoTransz(K::AbstractArray{T, 3}, d::NTuple{3, T}) where {T<:Number}
    Kreci = T(1) ./ K
    return T(2) .* prod(d) / d[3]^2 ./ (Kreci[:,:,1:end-1] + Kreci[:,:,2:end])
end
