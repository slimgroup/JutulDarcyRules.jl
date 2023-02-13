export setup_well_model

function force(M::jutulModel{D, T}, w::jutulForce{D, T}, tstep::Vector{T};
    ρCO2::T=T(ρCO2), ρH2O::T=T(ρH2O), g::T=T(10.0)) where {D, T}

    ## set up well information
    cell_loc = [Int.(round.(w.loc[i] ./ M.d)) for i = 1:length(w.loc)]
    Is = [setup_well(CartesianMesh(M), M.K, [cell_loc[i]], name = w.name[i]) for i = 1:length(w.loc)]
    ctrls = [w.name[i]==:Injector ? InjectorControl(TotalRateTarget(w.irate), [1.0, 0.0], density = ρCO2) : ProducerControl(BottomHolePressureTarget(2.0 * ρH2O * g * w.loc[i][3])) for i = 1:length(w.loc)]
    controls = Dict()
    for i = 1:length(w.loc)
        controls[w.name[i]] = ctrls[i]
    end
    return Is, controls
end

function setup_well_model(M::jutulModel{D, T}, f::jutulForce{D, T}, tstep::Vector{T};
    visCO2::T=T(visCO2), visH2O::T=T(visH2O), ρCO2::T=T(ρCO2), ρH2O::T=T(ρH2O), g::T=T(10.0)) where {D, T}

    ### set up well controls
    Is, controls = force(M, f, tstep; ρCO2=ρCO2, ρH2O=ρH2O, g=g)    

    ### set up model, parameters
    sys = ImmiscibleSystem((VaporPhase(), AqueousPhase()), reference_densities = [ρH2O, ρCO2])
    domain = discretized_domain_tpfv_flow(tpfv_geometry(CartesianMesh(M)), porosity = M.ϕ, permeability = M.K)
    model, parameters = setup_reservoir_model(domain, sys, wells = Is)
    select_output_variables!(model.models.Reservoir, :all)
    ρ = ConstantCompressibilityDensities(p_ref = 150*bar, density_ref = [ρCO2, ρH2O], compressibility = [1e-4/bar, 1e-6/bar])
    replace_variables!(model, PhaseMassDensities = ρ)
    replace_variables!(model, PhaseViscosities = vcat(visCO2 * ones(prod(M.n))', visH2O * ones(prod(M.n))'))
    replace_variables!(model, RelativePermeabilities = BrooksCoreyRelPerm(sys, [2.0, 2.0], [0.1, 0.1], 1.0))
    for x ∈ keys(model.models)
        Jutul.select_output_variables!(model.models[x], :all)
    end

    ### forces
    forces = setup_reservoir_forces(model, control = controls)

    ### initial state
    Z = repeat((1:M.n[end])*M.d[end], inner = prod(M.n[1:2]))
    p0 = ρH2O * g * Z # rho * g * h
    state0 = setup_reservoir_state(model, Pressure = 150*bar, Saturations = [0.0, 1.0])

    return model, parameters, state0, forces
end
