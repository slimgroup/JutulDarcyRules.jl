export jutulForce

struct jutulForce{D, T}
    irate::T
    name::Vector{Symbol}
    loc::Vector{NTuple{D, T}}
end

jutulForce(irate::T, loc::Vector{NTuple{D, T}}) = jutulForce(irate::T, vcat(:Injector, [:Producer for i = 1:length(loc)]), loc)

display(f::jutulForce{D, T}) where {D, T} =
    println("$(D)D jutulForce structure with $(length(f.loc)) injection/production wells and rate $(f.irate) m^3/s")

function force(M::jutulModel{D, T}, w::jutulForce{D, T}, tstep::Vector{T}; ρCO2::T=T(ρCO2)) where {D, T}
    model = model_(M)
    cell_loc = [Int.(round.(w.loc[i] ./ M.d)) for i = 1:length(w.loc)]
    Is = [setup_well(CartesianMesh(M), M.K, [cell_loc[i]], name = w.name[i]) for i = 1:length(w.loc)]
    ctrls = [w.name[i]==:Injector ? InjectorControl(TotalRateTarget(irate * ρCO2 * sum(tstep)), [1.0, 0.0], density = ρCO2) : ProducerControl(BottomHolePressureTarget(50*bar)) for i = 1:length(w.loc)]
    controls = Dict()
    for i = 1:length(w.loc)
        controls[w.name[i]] = ctrls[i]
    end
    return Is, controls
end

==(A::jutulForce{D, T}, B::jutulForce{D, T}) where {D,T} = (A.irate == B.irate && A.name == B.name && A.loc == B.loc)