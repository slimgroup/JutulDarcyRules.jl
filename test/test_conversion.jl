nx = 30
ny = 10
nz = 15
dims = (nx, ny, nz)
g = tpfv_geometry(CartesianMesh(dims, (rand(), rand(), rand())))
K1 = vcat(vec(rand(nx, ny, nz))', vec(rand(nx, ny, nz))', vec(rand(nx, ny, nz))')

## Test if same as JutulDarcy
@test compute_face_trans_(g, K1) == compute_face_trans(g, K1)
