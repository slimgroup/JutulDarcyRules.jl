@testset "Test jutulModel" begin
    
    @info "set up model"
    n = (rand(3:5), rand(3:5), rand(3:5))
    d = 10 .* (rand(), rand(), rand())
    K = 200 * rand() * md * ones(n)
    ϕ = rand()
    model = jutulModel(n, d, ϕ, K1to3(K))

    @test model.n == n
    @test model.d == d
    @test model.ϕ == ϕ
    @test model.K == K1to3(K)

    model1 = jutulModel(n, d, ϕ, K1to3(K))
    @test model1 == model
end