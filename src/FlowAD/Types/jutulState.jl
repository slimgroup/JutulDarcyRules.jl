export jutulState, jutulInitState, jutulStates, dict

abstract type jutulAllState{T} end

struct jutulState{T} <: jutulAllState{T}
    Saturations::Matrix{T}
    Pressure::Vector{T}
    TotalMasses::Matrix{T}
end

struct jutulStates{T} <: jutulAllState{T}
    states::Vector{jutulState{T}}
end

struct jutulInitState{T} <: jutulAllState{T}
    PhaseMassMobilities::Matrix{T}
    PhaseMassDensities::Matrix{T}
    RelativePermeabilities::Matrix{T}
    Saturations::Matrix{T}
    Pressure::Vector{T}
    TotalMasses::Matrix{T}
end

display(state::jutulAllState{T}) where T = println("$(typeof(state))")

jutulState(state::Dict{Symbol, T}) where T = jutulState(state[:Saturations], state[:Pressure], state[:TotalMasses])

jutulStates(states::Vector{Dict{Symbol, T}}) where T = jutulStates([jutulState(states[i]) for i = 1:length(states)])

jutulInitState(state::Dict{Symbol, T}) where T = jutulInitState(state[:PhaseMassMobilities],
    state[:PhaseMassDensities], state[:RelativePermeabilities], state[:Saturations],
    state[:Pressure], state[:TotalMasses])

function jutulInitState(M::jutulModel{D, T}; ρH2O::Number=T(1053.0), g::Number=T(10.0)) where {D, T}
    ## default state at time 0 with all water
    Z = repeat((1:M.n[end])*M.d[end], inner = prod(M.n[1:2]))
    p0 = ρH2O * g * Z # rho * g * h
    state0 = setup_state(model_(M), Pressure = p0, Saturations = [0.0, 1.0])
    return jutulInitState(state0)
end

dict(state::jutulAllState{T}) where T =
    Dict(fieldnames(typeof(state)) .=> getfield.(Ref(state), fieldnames(typeof(state))))

dict(state::jutulStates{T}) where T = [dict(state.states[i]) for i = 1:length(state.states)]
