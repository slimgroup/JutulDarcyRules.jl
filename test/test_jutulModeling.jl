model, model0, q, init_state, tstep = test_config();

## set up modeling operator
S = jutulModeling(model, tstep)

## simulation
x = log.(KtoTrans(CartesianMesh(model), model.K))
x0 = log.(KtoTrans(CartesianMesh(model0), model0.K))

states = S(x, q)

@testset "Test mass conservation" begin
    for i = 1:length(states.states)
        exist_co2 = sum(Saturations(states.states[i]) .* states.states[i].state[:Reservoir][:PhaseMassDensities][1,:] .* model.ϕ) * prod(model.d)
        inj_co2 = JutulDarcyAD.ρCO2 * q.irate * JutulDarcyAD.day * sum(tstep[1:i])
        @test isapprox(exist_co2, inj_co2) rtol=1e-3
    end
end