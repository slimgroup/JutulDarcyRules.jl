export jutulModeling

struct jutulModeling{D, T}
    model::jutulModel{D, T}
    state0::Dict
    tstep::Vector{T}
end

display(M::jutulModeling{D, T}) where {D, T} =
    println("$(M.D)D jutulModeling structure with:\n$(sum(M.tstep)) days in $(length(M.tstep)) time steps")

function (S::jutulModeling{D, T})(Transmissibilities::AbstractVector{T}, f::jutulForce{D, N};
    state0::jutulState=jutulState(S.model), visCO2::Number=T(1e-4), visH2O::Number=T(1e-3), ρCO2::Number=T(501.9), ρH2O::Number=T(1053.0)) where {D, T, N}

    model = model_(S.model)
    model.domain.grid.trans .= Transmissibilities
    forces = force(model, f; ρCO2=ρCO2)
    parameters = setup_parameters(model, PhaseViscosities = [visCO2, visH2O], density = [ρCO2, ρH2O]); # 0.1 and 1 cP
    states, rep = simulate(state0, model, day * S.tstep, parameters = parameters, forces = forces, info_level = 1, max_timestep_cuts = 1000)

end
