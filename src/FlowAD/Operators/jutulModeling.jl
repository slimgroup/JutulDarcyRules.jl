export jutulModeling

struct jutulModeling{D, T}
    model::jutulModel{D, T}
    tstep::Vector{T}
end

display(M::jutulModeling{D, T}) where {D, T} =
    println("$(D)D jutulModeling structure with $(sum(M.tstep)) days in $(length(M.tstep)) time steps")

function (S::jutulModeling{D, T})(Transmissibilities::AbstractVector{T}, f::jutulForce{D, N};
    state0::jutulInitState{T}=jutulInitState(S.model), visCO2::T=T(1e-4), visH2O::T=T(1e-3),
    ρCO2::T=T(501.9), ρH2O::T=T(1053.0), info_level::Int64=-1) where {D, T, N}

    forces = force(S.model, f; ρCO2=ρCO2)
    tstep = day * S.tstep
    model = model_(S.model)
    model.domain.grid.trans .= Transmissibilities
    parameters = setup_parameters(model, PhaseViscosities = [visCO2, visH2O], density = [ρCO2, ρH2O]); # 0.1 and 1 cP
    states, rep = simulate(dict(state0), model, tstep, parameters = parameters, forces = forces, info_level = info_level, max_timestep_cuts = 1000)
    return jutulStates(states)
end
