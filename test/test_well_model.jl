model, model0, q, q1, q2, init_state, init_state1, tstep = test_config();

@testset "Test parameters" begin
    visCO2 = 5e-5
    visH2O = 3e-3
    ρCO2 = 6e-7
    ρH2O = 7e-4
    jutul_model, parameters, state0, forces = setup_well_model(model, f, tstep; visCO2=visCO2, visH2O=visH2O, ρCO2=ρCO2, ρH2O=ρH2O)
    @test parameters[:Reservoir][:PhaseViscosities][1, :] .== visCO2
    @test parameters[:Reservoir][:PhaseViscosities][2, :] .== visH2O
end
