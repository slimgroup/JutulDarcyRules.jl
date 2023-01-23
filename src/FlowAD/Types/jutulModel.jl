export jutulModel

struct jutulModel{D, T}
    n::NTuple{D, Int64}
    d::NTuple{D, T}
    g::CartesianMesh{NTuple{D, Int64}, NTuple{T, T}, Vector{T}}
    ϕ::Union{Vector{T}, T}
    K::Union{Vector{T}, T}
end

function jutulModel(n::NTuple{D, Int64}, d::NTuple{D, T}, ϕ::Union{Vector{T}, T}, K::Union{Vector{T}, T})
    g = CartesianMesh(n, d .* n)
    return jutulModel(n, d, g, ϕ, K)
end

display(M::jutulModel{D, T}) where {D, T} =
    println("$(M.D)D jutulModel with size $(M.n) and grid spacing $(M.d)")

function model_(M::jutulModel{D, T}) where {D, T}
    g = CartesianMesh(M.n, M.d .* M.n)
    G = discretized_domain_tpfv_flow(tpfv_geometry(g), porosity = M.ϕ, permeability = M.K)
    model = SimulationModel(G, sys)
    replace_variables!(model, RelativePermeabilities = kr)
    return model
end