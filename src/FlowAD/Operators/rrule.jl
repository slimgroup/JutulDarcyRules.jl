function rrule(S::jutulModeling{D, T}, LogTransmissibilities::AbstractVector{T}, f::jutulForce{D, N};
    state0=nothing, visCO2::T=T(visCO2), visH2O::T=T(visH2O),
    ρCO2::T=T(ρCO2), ρH2O::T=T(ρH2O), info_level::Int64=-1) where {D, T, N}
    
    Transmissibilities = exp.(LogTransmissibilities)

    ### set up simulation time
    tstep = day * S.tstep

    ### set up simulation configurations
    model, parameters, state0_, forces = setup_well_model(S.model, f, tstep; visCO2=visCO2, visH2O=visH2O, ρCO2=ρCO2, ρH2O=ρH2O)
    model.models.Reservoir.domain.grid.trans .= Transmissibilities
    parameters[:Reservoir][:Transmissibilities] = Transmissibilities

    if isnothing(state0)
        state0 = state0_
    else
        state0 = dict(state0)
    end

    ### simulation
    sim, config = setup_reservoir_simulator(model, state0, parameters);
    states, report = simulate!(sim, tstep, forces = forces, config = config, max_timestep_cuts = 1000, info_level=info_level);
    output = jutulStates(states)
    
    ### optimization framework
    cfg = optimization_config(model, parameters, Dict(:Reservoir => [:FluidVolume, :Transmissibilities], :Injector => [:FluidVolume]))
    cfg[:Reservoir][:Transmissibilities][:scaler] = :log

    function pullback(dy)
        states_ref_ = output(vec(output)-dy)
        check_valid_state(states_ref_)
        states_ref = dict(states_ref_)
        mass_mismatch = (m, state, dt, step_no, forces) -> loss_per_step(m, state, dt, step_no, forces, states_ref)
        F_o, dF_o, F_and_dF, x0, lims, data = setup_parameter_optimization(
            model, state0, parameters, tstep, forces, mass_mismatch, cfg, param_obj = true, print = info_level, config = config, use_sparsity = false);
        g = similar(x0);
        misfit = F_and_dF(F_o, g, x0);
        return NoTangent(), g[1:length(LogTransmissibilities)], NoTangent()
    end
    return output, pullback
end

function loss_per_step(m, state, dt, step_no, forces, states_ref)
    state_ref = states_ref[step_no]
    fld = :Saturations
    fld2 = :Pressure
    val = state[:Reservoir][fld]
    val2 = state[:Reservoir][fld2]
    ref = state_ref[:Reservoir][fld]
    ref2 = state_ref[:Reservoir][fld2]
    return inner_mismatch(val, ref, val2, ref2)
end

function inner_mismatch(val, ref, val2, ref2)
    mismatch_s = zero(eltype(val))
    for i in axes(val, 2)
        mismatch_s += (val[1,i] - ref[1,i])^2
    end
    mismatch_p = zero(eltype(val2))
    for i in eachindex(val2)
        mismatch_p += (val2[i] - ref2[i])^2
    end
    return eltype(val)(0.5) * mismatch_s + eltype(val2)(0.5) * mismatch_p
end
