nx = 30
ny = 10
nz = 15
dims = (nx, ny, nz)
g = tpfv_geometry(CartesianMesh(dims, (rand(), rand(), rand())))
K = vcat(vec(rand(nx, ny, nz))', vec(rand(nx, ny, nz))', vec(rand(nx, ny, nz))')
@test compute_face_trans_(g, K) == compute_face_trans(g, K)
