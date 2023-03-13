model, model0, q, q1, q2, init_state, init_state1, tstep = test_config();

@testset "Test parameters" begin
    visCO2 = 5e-5
    visH2O = 3e-3
    ρCO2 = 6e-7
    ρH2O = 7e-4
    jutul_model, parameters, state0, forces = JutulDarcyRules.setup_well_model(model, q, tstep; visCO2=visCO2, visH2O=visH2O, ρCO2=ρCO2, ρH2O=ρH2O)
    @test all(parameters[:Reservoir][:PhaseViscosities][1, :] .== visCO2)
    @test all(parameters[:Reservoir][:PhaseViscosities][2, :] .== visH2O)
    @test jutul_model.models.Reservoir.system.rho_ref[1] == ρCO2
    @test jutul_model.models.Reservoir.system.rho_ref[2] == ρH2O

    @test jutul_model.models.Reservoir.secondary_variables[:PhaseMassDensities].reference_densities[1] == ρCO2
    @test jutul_model.models.Reservoir.secondary_variables[:PhaseMassDensities].reference_densities[2] == ρH2O
end
