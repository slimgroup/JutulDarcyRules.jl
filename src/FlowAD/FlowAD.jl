export jutulModeling

struct jutulModeling{T}
    model::SimulationModel
    parameters::Dict
    state0::Dict
    tstep::Vector{T}
    forces::NamedTuple
end

display(S::jutulModeling{T}) where T = println("jutulModeling structure with:\n$(sum(S.tstep)) days in $(length(S.tstep)) time steps")

function (S::jutulModeling{T})(K::AbstractArray{T}) where T
    states, rep = simulate(S.state0, S.model, S.tstep, parameters = S.parameters, forces = S.forces, info_level = 1, max_timestep_cuts = 1000)
end
