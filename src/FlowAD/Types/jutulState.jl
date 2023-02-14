export jutulState, jutulStates, dict, Saturations, Pressure

abstract type jutulAllState{T} <: DenseVector{T} end

struct jutulState{T} <: jutulAllState{T}
    state::Dict
end

struct jutulStates{T} <: jutulAllState{T}
    states::Vector{jutulState{T}}
end

display(state::jutulAllState{T}) where T = println("$(typeof(state))")

jutulState(state::Dict) = jutulState{eltype(state[:Reservoir][:Saturations])}(state)
jutulStates(states::Vector{Dict{Symbol, T}}) where T = jutulStates{eltype(states[1][:Reservoir][:Saturations])}([jutulState(states[i]) for i = 1:length(states)])

Saturations(state::jutulState) = state.state[:Reservoir][:Saturations][1,:]
Pressure(state::jutulState) = state.state[:Reservoir][:Pressure]

Saturations(state::jutulStates) = vcat([Saturations(state.states[i]) for i = 1:get_nt(state)]...)
Pressure(state::jutulStates) = vcat([Pressure(state.states[i]) for i = 1:get_nt(state)]...)

get_nt(state::jutulStates) = length(state.states)
get_nn(state::jutulState) = length(Saturations(state))
get_nn(state::jutulStates) = get_nn(state.states[1])

###### turn jutulStates to state dictionary

dict(state::jutulState) = state.state
dict(state::jutulStates) = [dict(state.states[i]) for i = 1:get_nt(state)]
dict(state::Dict) = state

###### AbstractVector

length(state::jutulState) = length(Saturations(state)) + length(Pressure(state))
length(state::jutulStates) = sum([length(state.states[i]) for i = 1:get_nt(state)])
size(state::jutulAllState) = (length(state),)

vec(state::jutulAllState) = vcat(Saturations(state), Pressure(state))

IndexStyle(::jutulAllState) = IndexLinear()
function getindex(state::jutulState, i::Int)
    if i <= get_nn(state)
        ## getindex for saturation
        return state.state[:Reservoir][:Saturations][1,i]
    else
        ## getindex for pressure
        return state.state[:Reservoir][:Pressure][i-get_nn(state)]
    end
end

function setindex!(state::jutulState{T}, v, i::Int) where T
    if i <= get_nn(state)
        ## setindex for saturation
        state.state[:Reservoir][:Saturations][1,i] = T(v)
        state.state[:Reservoir][:Saturations][2,i] = T(1) - T(v)
    else
        ## setindex for pressure
        state.state[:Reservoir][:Pressure][i-get_nn(state)] = T(v)
    end
end

function getindex(states::jutulStates, i::Int)
    idx_t = div(i-1, get_nn(states)) + 1
    idx_n = mod(i-1, get_nn(states)) + 1
    if idx_t <= get_nt(states)
        ## getindex for saturation
        return states.states[idx_t].state[:Reservoir][:Saturations][1,idx_n]
    else
        ## getindex for pressure
        return states.states[idx_t-get_nt(states)].state[:Reservoir][:Pressure][idx_n]
    end
end

function setindex!(states::jutulStates{T}, v, i::Int) where T
    idx_t = div(i-1, get_nn(states)) + 1
    idx_n = mod(i-1, get_nn(states)) + 1
    if idx_t <= get_nt(states)
        ## setindex for saturation
        states.states[idx_t].state[:Reservoir][:Saturations][1,idx_n] = T(v)
    else
        ## setindex for pressure
        states.states[idx_t-get_nt(states)].state[:Reservoir][:Pressure][idx_n] = T(v)
    end
end

getindex(state::jutulAllState, i...) = getindex(vec(state), i...)

firstindex(A::jutulAllState) = 1
lastindex(A::jutulAllState) = length(A)

norm(A::jutulAllState, order::Real=2) = norm(vec(A), order)

dot(A::jutulAllState, B::jutulAllState) = dot(vec(A), vec(B))
dot(A::jutulAllState, B::AbstractArray) = dot(vec(A), vec(B))
dot(A::AbstractArray, B::jutulAllState) = dot(vec(A), vec(B))

function (states::jutulState{T})(x::AbstractArray) where T
    @assert length(states) == length(x)
    states_ = deepcopy(states)
    states_ .= T.(x)
    states_.state[:Reservoir][:Saturations][2,:] = T(1) .- states_.state[:Reservoir][:Saturations][1,:]
    return states_
end

function (states::jutulStates{T})(x::AbstractArray) where T
    @assert length(states) == length(x)
    states_ = deepcopy(states)
    states_ .= T.(x)
    for i = 1:get_nt(states)
        states_.states[i].state[:Reservoir][:Saturations][2,:] = T(1) .- states_.states[i].state[:Reservoir][:Saturations][1,:]
    end
    return states_
end

for op in [:+, :-, :*, :/]
    @eval function $(op)(A::jutulAllState{T}, B::jutulAllState{T}) where T
        return A(broadcast($(op), vec(A), vec(B)))
    end
end

==(A::jutulAllState{T}, B::jutulAllState{T}) where {T} = vec(A) == vec(B)

function check_valid_state(states::jutulState{T}) where T
    @assert all(sum(states.state[:Reservoir][:Saturations], dims=1) .== 1)
end

function check_valid_state(states::jutulStates{T}) where T
    for i = 1:get_nt(states)
        check_valid_state(states.states[i])
    end
end