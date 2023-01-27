@testset "Test jutulForce" begin
    
    d = (1.0, 20.0, 3.0)
    inj_loc = (3, 1, 9) .* d
    prod_loc = (28, 1, 9) .* d
    irate = rand()
    q = jutulForce(irate, [inj_loc, prod_loc])
    q1 = jutulForce([irate, -irate], [inj_loc, prod_loc])
    @test q == q1
end