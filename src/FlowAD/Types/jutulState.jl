export jutulState, jutulInitState

abstract type jutulAllState{T} end

struct jutulState{T} <: jutulAllState{T}
    Saturations::Matrix{T}
    Pressure::Vector{T}
    TotalMasses::Matrix{T}
end

struct jutulInitState{T} <: jutulAllState{T}
    PhaseMassMobilities::Matrix{T}
    PhaseMassDensities::Matrix{T}
    RelativePermeabilities::Matrix{T}
    Saturations::Matrix{T}
    Pressure::Vector{T}
    TotalMasses::Matrix{T}
end

jutulState(state::Dict{Symbol, Any}) = jutulState(state[:Saturations], state[:Pressure], state[:TotalMasses])

jutulInitState(state::Dict{Symbol, Any}) = jutulInitState(state[:PhaseMassMobilities],
    state[:PhaseMassDensities], state[:RelativePermeabilities], state[:Saturations],
    state[:Pressure], state[:TotalMasses])

function jutulInitState(M::jutulModel{D, T}; ρH2O::Number=T(1053.0), g::Number=T(10.0)) where {D, T}
    ## default state at time 0 with all water
    Z = repeat((1:M.n[end])*M.d[end], inner = prod(M.n[1:2]))
    p0 = ρH2O * g * Z # rho * g * h
    state0 = setup_state(M, Pressure = p0, Saturations = [0.0, 1.0]);
    return jutulInitState(state0)
end

dict(S::jutulAllState{T}) where T = Dict(fieldnames(typeof(S)) .=> getfield.(Ref(S), fieldnames(typeof(S))))