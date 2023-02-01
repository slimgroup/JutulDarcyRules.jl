function rrule(S::jutulModeling{D, T}, LogTransmissibilities::AbstractVector{T}, f::jutulForce{D, N};
    state0::jutulState{T}=jutulState(S.model), visCO2::T=T(visCO2), visH2O::T=T(visH2O),
    ρCO2::T=T(ρCO2), ρH2O::T=T(ρH2O), info_level::Int64=-1) where {D, T, N}
    
    Transmissibilities = exp.(LogTransmissibilities)
    forces = force(S.model, f; ρCO2=ρCO2)
    tstep = day * S.tstep
    model = model_(S.model)
    model.domain.grid.trans .= Transmissibilities
    parameters = setup_parameters(model, PhaseViscosities = [visCO2, visH2O], density = [ρCO2, ρH2O]); # 0.1 and 1 cP
    states, _ = simulate(dict(state0), model, tstep, parameters = parameters, forces = forces, info_level = info_level, max_timestep_cuts = 1000)
    output = jutulStates(states)
    cfg = optimization_config(model, parameters, use_scaling = true, rel_min = 0.1, rel_max = 10)
    for (ki, vi) in cfg
        if ki in [:TwoPointGravityDifference, :PhaseViscosities]
            vi[:active] = false
        end
        if ki == :Transmissibilities
            vi[:scaler] = :log
        end
    end
    cfg[:Transmissibilities][:use_scaling] = false

    function pullback(dy)
        states_dy = output(dy)
        states_ref = dict(output-states_dy)
        function mass_mismatch(m, state, dt, step_no, forces)
            state_ref = states_ref[step_no]
            fld = :Saturations
            fld2 = :Pressure
            val = state[fld]
            val2 = state[fld2]
            ref = state_ref[fld]
            ref2 = state_ref[fld2]
            return 0.5 * sum((val[1,:] - ref[1,:]).^2) + 0.5 * sum((val2-ref2).^2)
        end
        Jutul.evaluate_objective(mass_mismatch, model, states_ref, tstep, forces)
        F_o, dF_o, F_and_dF, x0, lims, data = setup_parameter_optimization(model,
        dict(state0), parameters, tstep, forces, mass_mismatch, cfg, print = -1, param_obj = true);
        g = similar(x0);
        output = F_and_dF(F_o, g, x0);
        return NoTangent(), g[1:length(LogTransmissibilities)], NoTangent()
    end
    return output, pullback
end