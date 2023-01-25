export jutulModel, CartesianMesh

struct jutulModel{D, T}
    n::NTuple{D, Int64}
    d::NTuple{D, T}
    ϕ::Union{Vector{T}, T}
    K::Union{Matrix{T}, T}
end

display(M::jutulModel{D, T}) where {D, T} =
    println("$(D)D jutulModel with size $(M.n) and grid spacing $(M.d)")

CartesianMesh(M::jutulModel{D, T}) where {D, T} = CartesianMesh(M.n, M.d .* M.n)

function model_(M::jutulModel{D, T}) where {D, T}
    g = CartesianMesh(M.n, M.d .* M.n)
    G = discretized_domain_tpfv_flow(tpfv_geometry(g), porosity = M.ϕ, permeability = M.K)
    model = SimulationModel(G, sys)
    replace_variables!(model, RelativePermeabilities = kr)
    return model
end