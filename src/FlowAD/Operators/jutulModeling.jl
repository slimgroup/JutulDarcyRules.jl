export jutulModeling

struct jutulModeling{D, T}
    model::jutulModel{D, T}
    tstep::Vector{T}
end

display(M::jutulModeling{D, T}) where {D, T} =
    println("$(D)D jutulModeling structure with $(sum(M.tstep)) days in $(length(M.tstep)) time steps")

function (S::jutulModeling{D, T})(LogTransmissibilities::AbstractVector{T}, f::jutulForce{D, N};
    state0::jutulState{T}=jutulState(S.model), visCO2::T=T(visCO2), visH2O::T=T(visH2O),
    ρCO2::T=T(ρCO2), ρH2O::T=T(ρH2O), info_level::Int64=-1) where {D, T, N}

    Transmissibilities = exp.(LogTransmissibilities)

    sys = ImmiscibleSystem((VaporPhase(), AqueousPhase()), reference_densities = [ρH2O, ρCO2])

    ### set up model, parameters
    domain = discretized_domain_tpfv_flow(tpfv_geometry(CartesianMesh(S.model)), porosity = S.model.ϕ, permeability = S.model.K)
    model, parameters = setup_reservoir_model(domain, sys, wells = Is)
    select_output_variables!(model.models.Reservoir, :all)
    ρ = ConstantCompressibilityDensities(p_ref = 150*bar, density_ref = rhoS, compressibility = c)
    replace_variables!(model, PhaseMassDensities = ρ)
    replace_variables!(model, PhaseViscosities = [visCO2, visH2O])
    model.models.Reservoir.domain.grid.trans .= Transmissibilities

    ### set up well controls
    tstep = day * S.tstep
    Is, controls = force(S.model, f,tstep; ρCO2=ρCO2)    
    forces = setup_reservoir_forces(model, control = controls)
    
    ### simulation
    sim, config = setup_reservoir_simulator(model, state0, parameters);
    states, _ = simulate!(sim, tstep, forces = forces, config = config, info_level = info_level, max_timestep_cuts = 1000);

    return jutulStates(states)
end
