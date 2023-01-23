export jutulForce

struct jutulForce{D, T}
    irate::Vector{T}
    loc::Vector{NTuple{D, T}}
end

jutulForce(irate::T, loc::Vector{NTuple{D, T}}) = jutulForce(irate_ = [(-T(1))^(i+1)*irate for i = 1:length(loc)], loc)

display(f::jutulForce{D, T}) where {D, T}
    = println("$(f.D)D jutulForce structure with $(length(f.loc)) injection/production wells and rate $(f.irate)")

function force(M::jutulModel{D, T}, f::jutulForce{D, T}; ρCO2::Number=501.9) where {D, T}
    G = discretized_domain_tpfv_flow(tpfv_geometry(M.g), porosity = M.ϕ, permeability = M.K)
    model = SimulationModel(G, sys)
    cell_loc = [Int.(round.(f.loc[i] ./ M.d)) for i = 1:length(f.loc)]
    cell = [sum([(cell_loc[i][d]-1) * prod(M.n[1:d-1]) for d = length(cell_loc[i]):-1:1]) + 1 for i = 1:length(cell_loc)]
    src  = [SourceTerm(cell[i], f.irate[i] * ρCO2, fractional_flow = [T(f.irate[i] > 0), T(1)-T(f.irate[i] > 0)]) for i = 1:length(f.loc)]
    return setup_forces(model, sources = src)
end
