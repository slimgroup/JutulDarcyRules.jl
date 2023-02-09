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

function model_(M::jutulModel{D, T}; ρCO2::T=T(ρCO2), ρH2O::T=T(ρH2O)) where {D, T}
    g = CartesianMesh(M.n, M.d .* M.n)
    G = discretized_domain_tpfv_flow(tpfv_geometry(g), porosity = M.ϕ, permeability = M.K)
    model = SimulationModel(G, sys, output_level = :all)
    p_ref = 1e7 # 100 bar
    c = [1e-4, 1e-6]./1e5
    rho_at_ref = [ρCO2, ρH2O]
    dens = ConstantCompressibilityDensities(p_ref = p_ref, density_ref = rho_at_ref, compressibility = c)
    replace_variables!(model, PhaseMassDensities = dens)
    replace_variables!(model, RelativePermeabilities = BrooksCoreyRelPerm(sys, [2.0, 2.0], [0.1, 0.1], 1.0))
    return model
end

==(A::jutulModel{D, T}, B::jutulModel{D, T}) where {D,T} = (A.n == B.n && A.d == B.d && A.ϕ == B.ϕ && A.K == B.K)