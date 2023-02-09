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

    if isnothing(state0)
        state0 = state0_
    else
        state0 = dict(state0)
    end

    ### simulation
    sim, config = setup_reservoir_simulator(model, state0, parameters);
    states, report = simulate!(sim, tstep, forces = forces, config = config, max_timestep_cuts = 1000, info_level=info_level);
    output = jutulStates(states)
    return output
end
