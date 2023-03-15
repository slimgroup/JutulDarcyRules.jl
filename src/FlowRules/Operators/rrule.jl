function rrule(S::jutulModeling{D, T}, LogTransmissibilities::AbstractVector{T}, ϕ::AbstractVector{T}, f::Union{jutulForce{D, N}, jutulVWell{D, N}};
    state0=nothing, visCO2::T=T(visCO2), visH2O::T=T(visH2O),
    ρCO2::T=T(ρCO2), ρH2O::T=T(ρH2O), info_level::Int64=-1) where {D, T, N}
   
    Transmissibilities = exp.(LogTransmissibilities)

    ### set up simulation time
    tstep = day * S.tstep

    ### set up simulation configurations
    model, parameters, state0_, forces = setup_well_model(S.model, f, tstep; visCO2=visCO2, visH2O=visH2O, ρCO2=ρCO2, ρH2O=ρH2O)
    model.models.Reservoir.domain.grid.trans .= Transmissibilities
    model.models.Reservoir.domain.grid.pore_volumes .= prod(S.model.d) .* ϕ
    parameters[:Reservoir][:Transmissibilities] = Transmissibilities
    parameters[:Reservoir][:FluidVolume] .= prod(S.model.d) .* ϕ

    isnothing(state0) || (state0_[:Reservoir] = get_Reservoir_state(state0))

    ### simulation
    sim, config = setup_reservoir_simulator(model, state0_, parameters);
    states, reports = simulate!(sim, tstep, forces = forces, config = config, max_timestep_cuts = 1000, info_level=info_level);
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
            states, reports, model, state0_, parameters, tstep, forces, mass_mismatch, cfg, param_obj = true, print = info_level, config = config, use_sparsity = false);
        g = dF_o(similar(x0), x0);
        return NoTangent(), g[1:length(LogTransmissibilities)], g[length(LogTransmissibilities)+1:length(LogTransmissibilities)+prod(S.model.n)] * prod(S.model.d), NoTangent()
    end
    return output, pullback
end

function rrule(S::jutulModeling{D, T}, LogTransmissibilities::AbstractVector{T}, ϕ::AbstractVector{T}, f::jutulSource{D, N};
    state0=nothing, visCO2::T=T(visCO2), visH2O::T=T(visH2O),
    ρCO2::T=T(ρCO2), ρH2O::T=T(ρH2O), info_level::Int64=-1) where {D, T, N}
    
    Transmissibilities = exp.(LogTransmissibilities)

    forces = source(S.model, f; ρCO2=ρCO2)

    ### set up simulation time
    tstep = day * S.tstep
    model = simple_model(S.model; ρCO2=ρCO2, ρH2O=ρH2O)
    model.domain.grid.trans .= Transmissibilities
    model.domain.grid.pore_volumes .= prod(S.model.d) .* ϕ
    parameters = setup_parameters(model, PhaseViscosities = [visCO2, visH2O]);
    state0_ = jutulSimpleState(S.model)
    isnothing(state0) || (state0_ = state0)
    states, reports = simulate(dict(state0_), model, tstep, parameters = parameters, forces = forces, info_level = info_level, max_timestep_cuts = 1000)
    output = jutulSimpleStates(states)
    cfg = optimization_config(model, parameters, use_scaling = false, rel_min = 0., rel_max = Inf)
    for (ki, vi) in cfg
        if ki in [:TwoPointGravityDifference, :PhaseViscosities]
            vi[:active] = false
        end
        if ki == :Transmissibilities
            vi[:scaler] = :log
        end
    end

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
        mass_mismatch = (m, state, dt, step_no, forces) -> loss_per_step_simple(m, state, dt, step_no, forces, states_ref)
        Jutul.evaluate_objective(mass_mismatch, model, states_ref, tstep, forces)
        F_o, dF_o, F_and_dF, x0, lims, data = setup_parameter_optimization(states, reports, model,
        dict(state0_), parameters, tstep, forces, mass_mismatch, cfg, print = -1, param_obj = true);
        g = dF_o(similar(x0), x0);
        return NoTangent(), g[1:length(LogTransmissibilities)], g[length(LogTransmissibilities)+1:length(LogTransmissibilities)+prod(S.model.n)] * prod(S.model.d), NoTangent()
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

function loss_per_step_simple(m, state, dt, step_no, forces, states_ref)
    state_ref = states_ref[step_no]
    fld = :Saturations
    fld2 = :Pressure
    val = state[fld]
    val2 = state[fld2]
    ref = state_ref[fld]
    ref2 = state_ref[fld2]
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

function setup_parameter_optimization(precomputed_states, reports, model, state0, param, dt, forces, G, arg...; kwarg...)
    case = JutulCase(model, dt, forces, state0 = state0, parameters = param)
    return setup_parameter_optimization(precomputed_states, reports, case, G, arg...; kwarg...)
end

function setup_parameter_optimization(precomputed_states, reports, case::JutulCase, G, opt_cfg = optimization_config(case.model, case.parameters);
                                                            grad_type = :adjoint,
                                                            config = nothing,
                                                            print = 1,
                                                            copy_case = true,
                                                            param_obj = false,
                                                            use_sparsity = true,
                                                            kwarg...)
    if copy_case
        case = Jutul.duplicate(case)
    end
    # Pick active set of targets from the optimization config and construct a mapper
    (; model, state0, parameters) = case
    if print isa Bool
        if print
            print = 1
        else
            print = Inf
        end
    end
    verbose = print > 0 && isfinite(print)
    targets = optimization_targets(opt_cfg, model)
    if grad_type == :numeric
        @assert length(targets) == 1
        @assert model isa SimulationModel
    else
        @assert grad_type == :adjoint
    end
    mapper, = variable_mapper(model, :parameters, targets = targets, config = opt_cfg)
    lims = optimization_limits(opt_cfg, mapper, parameters, model)
    if verbose
        print_parameter_optimization_config(targets, opt_cfg, model)
    end
    x0 = vectorize_variables(model, parameters, mapper, config = opt_cfg)
    for k in eachindex(x0)
        low = lims[1][k]
        high = lims[2][k]
        @assert low <= x0[k] "Computed lower limit $low for parameter #$k was larger than provided x0[k]=$(x0[k])"
        @assert high >= x0[k] "Computer upper limit $hi for parameter #$k was smaller than provided x0[k]=$(x0[k])"
    end
    data = Dict()
    data[:n_objective] = 1
    data[:n_gradient] = 1
    data[:obj_hist] = zeros(0)

    sim = Simulator(case)
    if isnothing(config)
        config = simulator_config(sim; info_level = -1, kwarg...)
    elseif !verbose
        config[:info_level] = -1
        config[:end_report] = false
    end
    data[:sim] = sim
    data[:sim_config] = config

    if grad_type == :adjoint
        adj_storage = setup_adjoint_storage(model, state0 = state0,
                                                   parameters = parameters,
                                                   targets = targets,
                                                   use_sparsity = use_sparsity,
                                                   param_obj = param_obj)
        data[:adjoint_storage] = adj_storage
        grad_adj = zeros(adj_storage.n)
    else
        grad_adj = similar(x0)
    end
    data[:case] = case
    data[:grad_adj] = grad_adj
    data[:mapper] = mapper
    data[:G] = G
    data[:targets] = targets
    data[:mapper] = mapper
    data[:config] = opt_cfg
    data[:last_obj] = Inf
    data[:x_hash] = hash(x0)
    data[:states] = precomputed_states
    data[:reports] = reports
    F = x -> objective_opt!(x, data, print)
    dF = (dFdx, x) -> gradient_opt!(dFdx, x, data)
    F_and_dF = (F, dFdx, x) -> objective_and_gradient_opt!(F, dFdx, x, data, print)
    return (F! = F, dF! = dF, F_and_dF! = F_and_dF, x0 = x0, limits = lims, data = data)
end