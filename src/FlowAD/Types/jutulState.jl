export jutulState, jutulStates, dict

abstract type jutulAllState{T} <: DenseVector{T} end

struct jutulState{T} <: jutulAllState{T}
    PhaseMassMobilities::Matrix{T}
    PhaseMassDensities::Matrix{T}
    RelativePermeabilities::Matrix{T}
    Saturations::Vector{T}
    Pressure::Vector{T}
    TotalMasses::Matrix{T}
end

struct jutulStates{T} <: jutulAllState{T}
    states::Vector{jutulState{T}}
end

display(state::jutulAllState{T}) where T = println("$(typeof(state))")

jutulState(state::Dict{Symbol, T}) where T =
    jutulState(state[:PhaseMassMobilities], state[:PhaseMassDensities], state[:RelativePermeabilities], 
    state[:Saturations][1,:], state[:Pressure], state[:TotalMasses])

jutulStates(states::Vector{Dict{Symbol, T}}) where T = jutulStates([jutulState(states[i]) for i = 1:length(states)])

function jutulState(M::jutulModel{D, T}; ρCO2::T=T(ρCO2), ρH2O::T=T(ρH2O), g::T=T(10.0)) where {D, T}
    ## default state at time 0 with all water
    Z = repeat((1:M.n[end])*M.d[end], inner = prod(M.n[1:2]))
    p0 = ρH2O * g * Z # rho * g * h
    state0 = setup_state(model_(M; ρCO2=ρCO2, ρH2O=ρH2O), Pressure = p0, Saturations = [0.0, 1.0])
    return jutulState(state0)
end

get_nt(state::jutulStates) = length(state.states)
get_nn(state::jutulState) = length(state.Saturations)
get_nn(state::jutulStates) = get_nn(state.states[1])

###### turn jutulStates to state dictionary

function dict(state::jutulState{T}) where T
    dict_ = Dict(fieldnames(typeof(state)) .=> getfield.(Ref(state), fieldnames(typeof(state))))
    dict_[:Saturations] = hcat(dict_[:Saturations], T(1) .- dict_[:Saturations])'
    return dict_
end

dict(state::jutulStates{T}) where T = [dict(state.states[i]) for i = 1:get_nt(state)]

###### AbstractVector

length(state::jutulState) = length(state.Saturations) + length(state.Pressure)
length(state::jutulStates) = sum([length(state.states[i]) for i = 1:get_nt(state)])
size(state::jutulAllState) = (length(state),)

vec(state::jutulStates) = vcat(vcat([state.states[i].Saturations for i = 1:get_nt(state)]...), vcat([state.states[i].Pressure for i = 1:get_nt(state)]...))
vec(state::jutulState) = vcat(state.Saturations, state.Pressure)

IndexStyle(::jutulAllState) = IndexLinear()
function getindex(state::jutulState, i::Int)
    if i <= get_nn(state)
        ## getindex for saturation
        return state.Saturations[i]
    else
        ## getindex for pressure
        return state.Pressure[i-get_nn(state)]
    end
end

function getindex(state::jutulStates, i::Int)
    idx_t = div(i-1, get_nn(state)) + 1
    idx_n = mod(i-1, get_nn(state)) + 1
    if idx_t <= get_nt(state)
        ## getindex for saturation
        return state.states[idx_t].Saturations[idx_n]
    else
        ## getindex for pressure
        return state.states[idx_t-get_nt(state)].Pressure[idx_n]
    end
end

function setindex!(state::jutulState, v, i::Int)
    if i <= get_nn(state)
        ## setindex for saturation
        state.Saturations[i] = v
    else
        ## setindex for pressure
        state.Pressure[i-get_nn(state)] = v
    end
end

function setindex!(state::jutulStates, v, i::Int)
    idx_t = div(i-1, get_nn(state)) + 1
    idx_n = mod(i-1, get_nn(state)) + 1
    if idx_t <= get_nt(state)
        ## setindex for saturation
        state.states[idx_t].Saturations[idx_n] = v
    else
        ## setindex for pressure
        state.states[idx_t-get_nt(state)].Pressure[idx_n] = v
    end
end

firstindex(A::jutulAllState) = 1
lastindex(A::jutulAllState) = length(A)

norm(A::jutulAllState, order::Real=2) = norm(vec(A), order)

dot(A::jutulAllState, B::jutulAllState) = dot(vec(A), vec(B))
dot(A::jutulAllState, B::AbstractArray) = dot(vec(A), vec(B))
dot(A::AbstractArray, B::jutulAllState) = dot(vec(A), vec(B))

function (states::jutulAllState)(x::AbstractArray)
    @assert length(states) == length(x)
    states_ = deepcopy(states)
    states_ .= x
    return states_
end

for op in [:+, :-, :*, :/]
    @eval function $(op)(A::jutulAllState{T}, B::jutulAllState{T}) where T
        return A(broadcast($(op), vec(A), vec(B)))
    end
end

==(A::jutulAllState{T}, B::jutulAllState{T}) where {T} = vec(A) == vec(B)