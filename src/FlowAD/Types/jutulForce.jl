export jutulForce

struct jutulForce{D, T}
    irate::Vector{T}
    loc::Vector{NTuple{D, T}}
end

jutulForce(irate::T, loc::Vector{NTuple{D, T}}) where {D, T} = jutulForce([(-T(1))^(i+1)*irate for i = 1:length(loc)], loc)

display(f::jutulForce{D, T}) where {D, T} =
    println("$(D)D jutulForce structure with $(length(f.loc)) injection/production wells and rate $(f.irate) m^3/s")

function force(M::jutulModel{D, T}, f::jutulForce{D, T}; ρCO2::T=T(ρCO2)) where {D, T}
    model = model_(M)
    cell_loc = [Int.(round.(f.loc[i] ./ M.d)) for i = 1:length(f.loc)]
    cell = [sum([(cell_loc[i][d]-1) * prod(M.n[1:d-1]) for d = length(cell_loc[i]):-1:1]) + 1 for i = 1:length(cell_loc)]
    src  = [SourceTerm(cell[i], f.irate[i] * ρCO2, fractional_flow = [T(f.irate[i] > 0), T(1)-T(f.irate[i] > 0)]) for i = 1:length(f.loc)]
    return setup_forces(model, sources = src)
end

==(A::jutulForce{D, T}, B::jutulForce{D, T}) where {D,T} = (A.irate == B.irate && A.loc == B.loc)