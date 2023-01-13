export compute_face_trans_

### Can get rid of this clutter when https://github.com/FluxML/Zygote.jl/issues/1357 is addressed

function compute_face_trans_(g::JutulGeometry, K::AbstractArray{T, N}) where {T, N}
    T_hf = compute_half_face_trans_(g.cell_centroids, g.face_centroids, g.normals, g.areas, K, g.neighbors)
    return compute_face_trans_(T_hf, g.neighbors)
end

function compute_face_trans_(T_hf, N)
    faces = ChainRulesCore.@ignore_derivatives get_facepos(N)[1]
    nf = size(N, 2)
    T = [sum(1 ./T_hf[findall(faces .== i)]) for i = 1:nf]
    T = 1 ./T
    return T
end

function compute_half_face_trans_(cell_centroids, face_centroids, face_normals, face_areas, perm, N)
    ### adapted from https://github.com/sintefmath/Jutul.jl/blob/00df79186f9e26240c0ddbb6124c53ef4b27b7be/src/discretization/finite-volume.jl#L5
    nf = size(N, 2)
    nhf = 2*nf
    dim = size(cell_centroids, 1)

    T_hf = similar(cell_centroids, nhf)
    faces, facePos = ChainRulesCore.@ignore_derivatives get_facepos(N)
    nc = length(facePos)-1
    if isa(perm, AbstractFloat)
        perm = repeat([perm], 1, nc)
    else
        perm::AbstractVecOrMat
    end

    # Sanity check
    @assert(dim == 2 || dim == 3)
    # Check cell centroids
    @assert(size(cell_centroids, 1) == dim)
    @assert(size(cell_centroids, 2) == nc)
    # Check face centroids
    @assert(size(face_centroids, 1) == dim)
    @assert(size(face_centroids, 2) == nf)
    # Check normals
    @assert(size(face_normals, 1) == dim)
    @assert(size(face_normals, 2) == nf)
    # Check areas
    @assert(length(face_areas) == nf)
    # Check perm
    @assert(size(perm, 2) == nc)
    # Check N, just in case
    @assert(size(N, 2) == nf)
    return vcat([vcat([compute_half_trace_trans_(fpos, cell, faces, face_areas, perm, dim, face_centroids, cell_centroids, face_normals, N) for fpos = facePos[cell]:(facePos[cell+1]-1)]...) for cell = 1:nc]...)

end

function compute_half_trace_trans_(fpos, cell, faces, face_areas, perm, dim, face_centroids, cell_centroids, face_normals, N)
    face = faces[fpos]

    A = face_areas[face]
    K = expand_perm(perm[:, cell], dim)
    C = face_centroids[:, face] - cell_centroids[:, cell]
    Nn = face_normals[:, face]
    if N[2, face] == cell
        Nn = -Nn
    end
    return compute_half_face_trans(A, K, C, Nn)
end