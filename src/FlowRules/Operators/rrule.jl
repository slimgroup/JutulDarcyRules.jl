const _StateDict = OrderedDict{Symbol, Any}

function ordered_state_dicts(states, key::Symbol)
    filtered = _StateDict[]
    for s in states
        if s isa AbstractDict && haskey(s, key)
            push!(filtered, s isa _StateDict ? s : OrderedDict{Symbol, Any}(s))
        end
    end
    return filtered
end

_step_index(step_info) =
    step_info isa Integer ? step_info :
    (step_info isa AbstractDict && haskey(step_info, :step) ? step_info[:step] : step_info)

# Helper function to safely extract step index, handling potential offsets
function _safe_step_index(step_no, max_steps)
    step_ix = _step_index(step_no)
    if step_ix isa Integer
        # Handle potential 0-based indexing or offsets
        # Try both the direct index and index+1 (in case it's 0-based)
        # But first ensure it's within bounds
        if step_ix < 1
            # Might be 0-based, convert to 1-based
            step_ix = step_ix + 1
        elseif step_ix > max_steps
            # Might be 1-based but out of bounds, clamp it
            step_ix = max_steps
        end
        return max(1, min(step_ix, max_steps))
    else
        return max_steps
    end
end

function rrule(S::jutulModeling{D, T}, LogTransmissibilities::AbstractVector{T}, ϕ::AbstractVector{T}, f::Union{jutulForce{D, N}, jutulVWell{D, N}};
    state0=nothing, visCO2::T=T(visCO2), visH2O::T=T(visH2O),
    ρCO2::T=T(ρCO2), ρH2O::T=T(ρH2O), info_level::Int64=-1) where {D, T, N}
   
    Transmissibilities = exp.(LogTransmissibilities)

    ### set up simulation time
    tstep = day * S.tstep

    ### set up simulation configurations
    model, parameters, state0_, forces = setup_well_model(S.model, f, tstep; visCO2=visCO2, visH2O=visH2O, ρCO2=ρCO2, ρH2O=ρH2O)

    model.models.Reservoir.data_domain[:porosity] = ϕ
    parameters[:Reservoir][:Transmissibilities] = Transmissibilities
    parameters[:Reservoir][:FluidVolume] .= prod(S.model.d) .* ϕ

    isnothing(state0) || (state0_[:Reservoir] = get_Reservoir_state(state0))

    ### simulation
    sim, config = setup_reservoir_simulator(model, state0_, parameters);
    states, reports = simulate!(sim, tstep, forces = forces, config = config, max_timestep_cuts = 1000, info_level=info_level);
    output = jutulStates(states)
    reservoir_states = ordered_state_dicts(states, :Reservoir)
    
    ### optimization framework
    cfg = optimization_config(model, parameters, Dict(:Reservoir => [:FluidVolume, :Transmissibilities], :Injector => [:FluidVolume]))
    cfg[:Reservoir][:Transmissibilities][:scaler] = :log

    function pullback(dy)
        states_ref_ = output(vec(output)-dy)
        check_valid_state(states_ref_)
        states_ref = dict(states_ref_)
        # Critical fix: Ensure states_ref uses the same filtering as reservoir_states
        # The key insight: dict(jutulStates) returns states in the same order as jutulStates.states
        # which should match the order of reservoir_states (both filtered from the same source)
        # However, we need to ensure they are filtered identically
        states_ref_filtered = ordered_state_dicts(states_ref, :Reservoir)
        @assert length(states_ref_filtered) == length(reservoir_states) "states_ref_filtered length $(length(states_ref_filtered)) != reservoir_states length $(length(reservoir_states))"
        # Use loss_per_step with the filtered states_ref to ensure correct indexing
        # We've ensured states_ref_filtered and reservoir_states have the same length and order
        # step_no should correspond to the index in reservoir_states
        mass_mismatch = (m, state, dt, step_no, forces) -> loss_per_step(m, state, dt, step_no, forces, states_ref_filtered)
        F_o, dF_o, F_and_dF, x0, lims, data = setup_parameter_optimization(
            reservoir_states, reports, model, state0_, parameters, tstep, forces, mass_mismatch, cfg, param_obj = true, print = info_level, config = config, use_sparsity = false);
        g = dF_o(similar(x0), x0);
        n_faces = length(LogTransmissibilities)
        n_cells = prod(S.model.n)
        dLogTransmissibilities = g[1:n_faces]
        dϕ = g[n_faces + 1 : n_faces + n_cells] * prod(S.model.d)
        return NoTangent(), dLogTransmissibilities, dϕ, NoTangent()
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
    model.data_domain[:porosity] = ϕ

    parameters = setup_parameters(model, PhaseViscosities = [visCO2, visH2O]);
    parameters[:Transmissibilities] = Transmissibilities
    parameters[:FluidVolume] .= prod(S.model.d) .* ϕ

    state0_ = jutulSimpleState(S.model)
    isnothing(state0) || (state0_ = state0)
    states, reports = simulate(dict(state0_), model, tstep, parameters = parameters, forces = forces, info_level = info_level, max_timestep_cuts = 1000)
    output = jutulSimpleStates(states)
    simple_states = ordered_state_dicts(states, :Saturations)
    cfg = optimization_config(model, parameters, use_scaling = false, rel_min = 0., rel_max = nothing)
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
        # Critical fix: Ensure states_ref uses the same filtering as simple_states
        # Both should be filtered from the same source using ordered_state_dicts
        # This ensures they have the same order and structure
        states_ref_filtered = ordered_state_dicts(states_ref, :Saturations)
        @assert length(states_ref_filtered) == length(simple_states) "states_ref_filtered length $(length(states_ref_filtered)) != simple_states length $(length(simple_states))"
        # Convert to the format expected by loss_per_step_simple
        # Use the filtered states to ensure correct ordering
        states_ref = OrderedDict{Symbol, Any}[OrderedDict{Symbol, Any}(s) for s in states_ref_filtered]
        mass_mismatch = (m, state, dt, step_no, forces) -> loss_per_step_simple(m, state, dt, step_no, forces, states_ref)
        Jutul.evaluate_objective(mass_mismatch, model, states_ref, tstep, forces)
        F_o, dF_o, F_and_dF, x0, lims, data = setup_parameter_optimization(simple_states, reports, model,
        dict(state0_), parameters, tstep, forces, mass_mismatch, cfg, print = -1, param_obj = true);
        g = dF_o(similar(x0), x0);
        n_faces = length(LogTransmissibilities)
        n_cells = prod(S.model.n)
        dLogTransmissibilities = g[1:n_faces]
        dϕ = g[n_faces + 1 : n_faces + n_cells] * prod(S.model.d)
        return NoTangent(), dLogTransmissibilities, dϕ, NoTangent()
    end
    return output, pullback
end

function loss_per_step(m, state, dt, step_no, forces, states_ref, reservoir_states=nothing)
    # Use the safe step index function to handle potential offsets
    step_ix = _safe_step_index(step_no, length(states_ref))
    # Critical: Use the step index to get the corresponding reference state
    # step_ix should correspond to the index in reservoir_states
    # states_ref should be in the same order as reservoir_states
    # If reservoir_states is provided, we can verify the index is correct
    if reservoir_states !== nothing
        # Ensure step_ix is within bounds of both arrays
        step_ix = max(1, min(step_ix, min(length(states_ref), length(reservoir_states))))
    end
    state_ref = states_ref[step_ix]
    fld = :Saturations
    fld2 = :Pressure
    val = state[:Reservoir][fld]
    val2 = state[:Reservoir][fld2]
    ref = state_ref[:Reservoir][fld]
    ref2 = state_ref[:Reservoir][fld2]
    return inner_mismatch(val, ref, val2, ref2)
end

function loss_per_step_simple(m, state, dt, step_no, forces, states_ref)
    # Use the safe step index function to handle potential offsets
    step_ix = _safe_step_index(step_no, length(states_ref))
    state_ref = states_ref[step_ix]
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
