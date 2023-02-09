model, model0, q, init_state, tstep = test_config();

## set up modeling operator
S = jutulModeling(model0, tstep)

## simulation
x = log.(KtoTrans(CartesianMesh(model), model.K))
x0 = log.(KtoTrans(CartesianMesh(model0), model0.K))

states = S(x, q)

misfit(x0) = 0.5 * norm(S(x0, q) - states).^2
g = gradient(()->misfit(x0), Flux.params(x0))

dx = randn(MersenneTwister(2023), length(x0))
dx = dx/norm(dx) * norm(x0)/50.0

@testset "Taylor-series gradient test of jutulModeling" begin
    grad_test(misfit, x0, dx, g[x0])
end