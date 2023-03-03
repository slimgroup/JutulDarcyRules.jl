export jutulModeling

struct jutulModeling{D, T}
    model::jutulModel{D, T}
    tstep::Vector{T}
end

display(M::jutulModeling{D, T}) where {D, T} =
    println("$(D)D jutulModeling structure with $(sum(M.tstep)) days in $(length(M.tstep)) time steps")

function (S::jutulModeling{D, T})(LogTransmissibilities::AbstractVector{T}, f::jutulForce{D, N};
    state0=nothing, visCO2::T=T(visCO2), visH2O::T=T(visH2O),
    ρCO2::T=T(ρCO2), ρH2O::T=T(ρH2O), info_level::Int64=-1) where {D, T, N}

    Transmissibilities = exp.(LogTransmissibilities)

    ### set up simulation time
    tstep = day * S.tstep

    ### set up simulation configurations
    model, parameters, state0_, forces = setup_well_model(S.model, f, tstep; visCO2=visCO2, visH2O=visH2O, ρCO2=ρCO2, ρH2O=ρH2O)
    model.models.Reservoir.domain.grid.trans .= Transmissibilities
    parameters[:Reservoir][:Transmissibilities] = Transmissibilities

    isnothing(state0) || (state0_[:Reservoir] = get_Reservoir_state(state0))

    ### simulation
    sim, config = setup_reservoir_simulator(model, state0_, parameters);
    states, report = simulate!(sim, tstep, forces = forces, config = config, max_timestep_cuts = 1000, info_level=info_level);
    output = jutulStates(states)
    return output
end

function (S::jutulModeling{D, T})(LogTransmissibilities::AbstractVector{T}, f::jutulSource{D, N};
    state0::jutulSimpleState{T}=jutulSimpleState(S.model), visCO2::T=T(visCO2), visH2O::T=T(visH2O),
    ρCO2::T=T(ρCO2), ρH2O::T=T(ρH2O), info_level::Int64=-1) where {D, T, N}

    Transmissibilities = exp.(LogTransmissibilities)

    forces = source(S.model, f; ρCO2=ρCO2)

    ### set up simulation time
    tstep = day * S.tstep
    model = simple_model(S.model; ρCO2=ρCO2, ρH2O=ρH2O)
    model.domain.grid.trans .= Transmissibilities
    parameters = setup_parameters(model, PhaseViscosities = [visCO2, visH2O]);
    states, _ = simulate(dict(state0), model, tstep, parameters = parameters, forces = forces, info_level = info_level, max_timestep_cuts = 1000)
    return jutulSimpleStates(states)
end