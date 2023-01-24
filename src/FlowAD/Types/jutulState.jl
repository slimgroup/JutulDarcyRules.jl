export jutulState, jutulInitState, jutulStates, dict

abstract type jutulAllState{T} <: DenseVector{T} end

struct jutulState{T} <: jutulAllState{T}
    Saturations::Vector{T}
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
    Saturations::Vector{T}
    Pressure::Vector{T}
    TotalMasses::Matrix{T}
end

display(state::jutulAllState{T}) where T = println("$(typeof(state))")

jutulState(state::Dict{Symbol, T}) where T = jutulState(state[:Saturations][1,:], state[:Pressure], state[:TotalMasses])

jutulStates(states::Vector{Dict{Symbol, T}}) where T = jutulStates([jutulState(states[i]) for i = 1:length(states)])

jutulInitState(state::Dict{Symbol, T}) where T = jutulInitState(state[:PhaseMassMobilities],
    state[:PhaseMassDensities], state[:RelativePermeabilities], state[:Saturations][1,:],
    state[:Pressure], state[:TotalMasses])

function jutulInitState(M::jutulModel{D, T}; ρH2O::T=T(1053.0), g::T=T(10.0)) where {D, T}
    ## default state at time 0 with all water
    Z = repeat((1:M.n[end])*M.d[end], inner = prod(M.n[1:2]))
    p0 = ρH2O * g * Z # rho * g * h
    state0 = setup_state(model_(M), Pressure = p0, Saturations = [0.0, 1.0])
    return jutulInitState(state0)
end

get_nt(state::jutulStates) = length(state.states)

###### turn jutulStates to state dictionary

function dict(state::jutulAllState{T}) where T
    dict_ = Dict(fieldnames(typeof(state)) .=> getfield.(Ref(state), fieldnames(typeof(state))))
    dict_[:Saturations] = hcat(dict_[:Saturations], T(1) .- dict_[:Saturations])'
    return dict_
end

dict(state::jutulStates{T}) where T = [dict(state.states[i]) for i = 1:get_nt(state)]

###### AbstractVector

length(state::jutulAllState) = length(state.Saturations) + length(state.Pressure)
length(state::jutulStates) = sum([length(state.states[i]) for i = 1:get_nt(state)])
size(state::jutulAllState) = (length(state),)

vec(state::jutulStates) = vcat([vcat(state.states[i].Saturations, state.states[i].Pressure) for i = 1:get_nt(state)]...)
vec(state::jutulAllState) = vcat(state.Saturations, state.Pressure)

IndexStyle(state::jutulAllState) = IndexLinear()
getindex(state::jutulAllState, i::Int) = vec(state)[i]
setindex!(state::jutulAllState{T}, v::T, i::Int) where T = begin vcat(state.Saturations, state.Pressure)[i] = v end
setindex!(state::jutulStates{T}, v::T, i::Int) where T = begin vcat([vcat(state.states[i].Saturations, state.states[i].Pressure) for i = 1:get_nt(state)]...)[i] = v end

firstindex(A::jutulAllState) = 1
lastindex(A::jutulAllState) = length(A)

norm(A::jutulAllState, order::Real=2) = norm(vec(A), order)

dot(A::jutulAllState, B::jutulAllState) = dot(vec(A), vec(B))
dot(A::jutulAllState, B::AbstractArray) = dot(vec(A), vec(B))
dot(A::AbstractArray, B::jutulAllState) = dot(vec(A), vec(B))
