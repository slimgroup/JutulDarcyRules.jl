nx = 30
ny = 10
nz = 15
dims = (nx, ny, nz)
g = tpfv_geometry(CartesianMesh(dims, (rand(), rand(), rand())))
K1 = vcat(vec(rand(nx, ny, nz))', vec(rand(nx, ny, nz))', vec(rand(nx, ny, nz))')

@testset "Test conversion" begin
    @info "compute transmissibility from permeability"
    @test compute_face_trans_(g, K1) == compute_face_trans(g, K1)
end