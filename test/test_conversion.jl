nx = 30
ny = 10
nz = 15
dims = (nx, ny, nz)
g1 = CartesianMesh(dims, (5.0, 8.0, 10.0) .* dims)
g = tpfv_geometry(g1)
K1 = vcat(vec(rand(nx, ny, nz))', vec(rand(nx, ny, nz))', vec(rand(nx, ny, nz))')

@testset "Test conversion" begin
    @info "compute transmissibility from permeability"
    @test isapprox(KtoTrans(g1, K1), compute_face_trans(g, K1))

    trans = KtoTrans(g1, K1)
    @info "compute permeability from transmissibility"
    Kinv = TransToK(g1, trans)
    @test isapprox(KtoTrans(g1, Kinv), trans) rtol=5e-3

end
