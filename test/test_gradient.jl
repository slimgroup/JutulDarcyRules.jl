model, model0, q, tstep = test_config()

## set up modeling operator
S = jutulModeling(model0, tstep)

## simulation
x = log.(KtoTrans(CartesianMesh(model), model.K))
x0 = log.(KtoTrans(CartesianMesh(model0), model0.K))

states = S(x, q)

dx = randn(length(x0))
dx = dx/norm(dx) * norm(x0)/100

misfit(x0) = 0.5 * norm(S(x0, q) - states).^2
g = gradient(()->misfit(x0), Flux.params(x0))

@info "test AD of jutulModeling"
@testset "Test AD of jutulModeling" begin
    grad_test(misfit, x0, dx, g[x0])
end