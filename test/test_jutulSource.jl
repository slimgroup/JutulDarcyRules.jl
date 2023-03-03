@testset "Test jutulSource" begin
    
    d = (1.0, 20.0, 3.0)
    inj_loc = (3, 1, 9) .* d
    prod_loc = (28, 1, 9) .* d
    irate = rand()
    q = jutulForce(irate, [inj_loc, prod_loc])
    q1 = jutulForce(irate, [inj_loc])

    @test q != q1
    @test q.irate == q1.irate
    @test q.loc[1] == q1.loc[1]
    @test q.name[1] == q1.name[1]
    @test q.name[2] == :Producer
end