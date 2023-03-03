model, model0, q, q1, init_state, init_state1, tstep = test_config();

## set up modeling operator
S = jutulModeling(model0, tstep)

## simulation
x = log.(KtoTrans(CartesianMesh(model), model.K))
x0 = log.(KtoTrans(CartesianMesh(model0), model0.K))

states = S(x, q)

misfit(x0) = 0.5 * norm(S(x0, q) - states).^2
g = gradient(()->misfit(x0), Flux.params(x0))

dx = randn(MersenneTwister(2023), length(x0))
dx = dx/norm(dx) * norm(x0)/5.0

@testset "Taylor-series gradient test of jutulModeling with wells" begin
    grad_test(misfit, x0, dx, g[x0])
end

states1 = S(x, q1)

misfit1(x0) = 0.5 * norm(S(x0, q1) - states1).^2
g1 = gradient(()->misfit1(x0), Flux.params(x0))

@testset "Taylor-series gradient test of simple jutulModeling" begin
    grad_test(misfit1, x0, dx/1.5, g1[x0])
end
