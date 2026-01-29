const _StateDict = OrderedDict{Symbol, Any}

"""
    ordered_state_dicts(states, key::Symbol)

Filter states to extract only those containing the specified key, maintaining order.
"""
function ordered_state_dicts(states, key::Symbol)
    filtered = _StateDict[]
    for s in states
        if s isa AbstractDict && haskey(s, key)
            push!(filtered, s isa _StateDict ? s : OrderedDict{Symbol, Any}(s))
        end
    end
    return filtered
end

"""
    _step_index(step_info)

Extract step index from step_info, handling both integer and dictionary formats.
"""
_step_index(step_info) =
    step_info isa Integer ? step_info :
    (step_info isa AbstractDict && haskey(step_info, :step) ? step_info[:step] : step_info)

"""
    _safe_step_index(step_no, max_steps)

Safely extract step index and clamp to valid range [1, max_steps].
"""
function _safe_step_index(step_no, max_steps)
    step_ix = _step_index(step_no)
    if step_ix isa Integer
        return clamp(step_ix, 1, max_steps)
    else
        return max_steps
    end
end

function rrule(S::jutulModeling{D, T}, LogTransmissibilities::AbstractVector{T}, ϕ::AbstractVector{T}, f::Union{jutulForce{D, N}, jutulVWell{D, N}};
    state0=nothing, visCO2::T=T(visCO2), visH2O::T=T(visH2O),
    ρCO2::T=T(ρCO2), ρH2O::T=T(ρH2O), info_level::Int64=-1) where {D, T, N}
   
    Transmissibilities = exp.(LogTransmissibilities)
    tstep = day * S.tstep

    model, parameters, state0_, forces = setup_well_model(S.model, f, tstep; visCO2=visCO2, visH2O=visH2O, ρCO2=ρCO2, ρH2O=ρH2O)

    model.models.Reservoir.data_domain[:porosity] = ϕ
    parameters[:Reservoir][:Transmissibilities] = Transmissibilities
    parameters[:Reservoir][:FluidVolume] .= prod(S.model.d) .* ϕ

    isnothing(state0) || (state0_[:Reservoir] = get_Reservoir_state(state0))

    sim, config = setup_reservoir_simulator(model, state0_, parameters);
    states, reports = simulate!(sim, tstep, forces = forces, config = config, max_timestep_cuts = 1000, info_level=info_level);
    output = jutulStates(states)
    reservoir_states = ordered_state_dicts(states, :Reservoir)
    
    # Configure optimization: only optimize Reservoir parameters
    cfg = optimization_config(model, parameters, Dict(:Reservoir => [:FluidVolume, :Transmissibilities]))
    cfg[:Reservoir][:Transmissibilities][:scaler] = :log

    function pullback(dy)
        # dy is dJ/d(output), where output = jutulStates from forward simulation
        # We need to compute dJ/d(LogTransmissibilities) and dJ/d(ϕ)
        # 
        # For VJP: dJ/d(params) = (dS/d(params))^T * dy
        # We construct a LINEAR objective: G(state) = <dy, state>
        # Then dG/d(params) = (dS/d(params))^T * dy, which is exactly the VJP!
        
        # Convert dy to jutulStates format to extract saturation and pressure components
        states_dy = output(dy)
        
        # Create weighted objective: G = Σ dy_i * state_i (linear inner product)
        function weighted_objective(m, state, dt, step_info, forces)
            step_ix = _safe_step_index(step_info, length(states_dy.states))
            dy_state = states_dy.states[step_ix]
            
            sat = state[:Reservoir][:Saturations]
            pres = state[:Reservoir][:Pressure]
            dy_sat = dy_state.state[:Reservoir][:Saturations]
            dy_pres = dy_state.state[:Reservoir][:Pressure]
            
            obj = zero(eltype(sat))
            for i in axes(sat, 2)
                obj += dy_sat[1, i] * sat[1, i]
            end
            for i in eachindex(pres)
                obj += dy_pres[i] * pres[i]
            end
            return obj
        end
        
        # Use solve_adjoint_sensitivities for Transmissibilities
        grad_dict = Jutul.solve_adjoint_sensitivities(
            model, states, reports, weighted_objective,
            forces = forces, state0 = state0_, parameters = parameters,
            raw_output = false
        )
        
        # Extract Transmissibilities gradient: dG/dT
        # Chain rule: dG/d(log(T)) = dG/dT * T
        dT = grad_dict[:Reservoir][:Transmissibilities]
        dLogTransmissibilities = dT .* Transmissibilities
        
        # For porosity, combine adjoint gradient (FluidVolume) with correction for 
        # direct porosity effect in data_domain
        # FluidVolume = cell_volume * ϕ, so dG/dϕ_FluidVolume = dG/dFluidVolume * cell_volume
        cell_volume = prod(S.model.d)
        dFluidVolume = grad_dict[:Reservoir][:FluidVolume]
        dϕ_from_FV = dFluidVolume .* cell_volume
        
        # Add contribution from direct porosity effect (through data_domain[:porosity])
        # This is computed via numerical differentiation on a few representative cells
        # and used to estimate a correction factor
        timesteps_sim = Jutul.report_timesteps(reports)
        n_cells = length(ϕ)
        dϕ = zeros(T, n_cells)
        epsilon = T(1e-4)
        
        # Compute base objective
        G_base = zero(T)
        for (step_ix, (s, dt)) in enumerate(zip(states, timesteps_sim))
            step_info = Dict{Symbol, Any}(:step => step_ix, :ministep => 1)
            G_base += weighted_objective(model, s, dt, step_info, forces)
        end
        
        # Compute gradient for each cell using central difference for better accuracy
        for i in 1:n_cells
            if ϕ[i] >= 1.0
                dϕ[i] = zero(T)
                continue
            end
            
            # Forward perturbation: ϕ[i] + epsilon
            ϕ_plus = copy(ϕ)
            ϕ_plus[i] += epsilon
            
            parameters_plus = deepcopy(parameters)
            parameters_plus[:Reservoir][:FluidVolume] = cell_volume .* ϕ_plus
            
            model_plus = deepcopy(model)
            model_plus.models.Reservoir.data_domain[:porosity] = ϕ_plus
            
            sim_plus, config_plus = setup_reservoir_simulator(model_plus, state0_, parameters_plus);
            states_plus, _ = simulate!(sim_plus, tstep, forces = forces, config = config_plus, max_timestep_cuts = 1000, info_level=-1);
            
            G_plus = zero(T)
            for (step_ix, (s, dt)) in enumerate(zip(states_plus, timesteps_sim))
                step_info = Dict{Symbol, Any}(:step => step_ix, :ministep => 1)
                G_plus += weighted_objective(model_plus, s, dt, step_info, forces)
            end
            
            # Backward perturbation: ϕ[i] - epsilon
            ϕ_minus = copy(ϕ)
            ϕ_minus[i] -= epsilon
            
            parameters_minus = deepcopy(parameters)
            parameters_minus[:Reservoir][:FluidVolume] = cell_volume .* ϕ_minus
            
            model_minus = deepcopy(model)
            model_minus.models.Reservoir.data_domain[:porosity] = ϕ_minus
            
            sim_minus, config_minus = setup_reservoir_simulator(model_minus, state0_, parameters_minus);
            states_minus, _ = simulate!(sim_minus, tstep, forces = forces, config = config_minus, max_timestep_cuts = 1000, info_level=-1);
            
            G_minus = zero(T)
            for (step_ix, (s, dt)) in enumerate(zip(states_minus, timesteps_sim))
                step_info = Dict{Symbol, Any}(:step => step_ix, :ministep => 1)
                G_minus += weighted_objective(model_minus, s, dt, step_info, forces)
            end
            
            # Central difference: more accurate, O(epsilon^2) error
            dϕ[i] = (G_plus - G_minus) / (T(2) * epsilon)
        end
        
        if info_level >= 0
            println("DEBUG: dLogT[1:5] = ", dLogTransmissibilities[1:min(5, length(dLogTransmissibilities))])
            println("DEBUG: dϕ[1:5] = ", dϕ[1:min(5, length(dϕ))])
        end
        
        return NoTangent(), dLogTransmissibilities, dϕ, NoTangent()
    end
    return output, pullback
end

function rrule(S::jutulModeling{D, T}, LogTransmissibilities::AbstractVector{T}, ϕ::AbstractVector{T}, f::jutulSource{D, N};
    state0=nothing, visCO2::T=T(visCO2), visH2O::T=T(visH2O),
    ρCO2::T=T(ρCO2), ρH2O::T=T(ρH2O), info_level::Int64=-1) where {D, T, N}
    
    Transmissibilities = exp.(LogTransmissibilities)
    forces = source(S.model, f; ρCO2=ρCO2)
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
        # dy is dJ/d(output), where output = jutulSimpleStates from forward simulation
        # We need to compute dJ/d(LogTransmissibilities) and dJ/d(ϕ)
        # 
        # For VJP: dJ/d(params) = (dS/d(params))^T * dy
        # We construct a LINEAR objective: G(state) = <dy, state>
        # Then dG/d(params) = (dS/d(params))^T * dy, which is exactly the VJP!
        
        # Convert dy to jutulSimpleStates format
        states_dy = output(dy)
        
        # Create weighted objective: G = Σ dy_i * state_i (linear inner product)
        function weighted_objective(m, state, dt, step_info, forces)
            step_ix = _safe_step_index(step_info, length(states_dy.states))
            dy_state = states_dy.states[step_ix]
            
            sat = state[:Saturations]
            pres = state[:Pressure]
            dy_sat = dy_state.state[:Saturations]
            dy_pres = dy_state.state[:Pressure]
            
            obj = zero(eltype(sat))
            for i in axes(sat, 2)
                obj += dy_sat[1, i] * sat[1, i]
            end
            for i in eachindex(pres)
                obj += dy_pres[i] * pres[i]
            end
            return obj
        end
        
        # Update parameters to match current values used in forward simulation
        parameters[:Transmissibilities] = Transmissibilities
        parameters[:FluidVolume] .= prod(S.model.d) .* ϕ
        model.data_domain[:porosity] = ϕ
        
        # Use solve_adjoint_sensitivities for Transmissibilities (it works correctly)
        grad_dict = Jutul.solve_adjoint_sensitivities(
            model, states, reports, weighted_objective,
            forces = forces, state0 = dict(state0_), parameters = parameters,
            raw_output = false
        )
        
        # Transmissibilities gradient: dG/dT
        # Chain rule: dG/d(log(T)) = dG/dT * T
        dT = grad_dict[:Transmissibilities]
        dLogTransmissibilities = dT .* Transmissibilities
        
        # For porosity, compute gradient via numerical differentiation to capture
        # all effects (both through FluidVolume and direct data_domain[:porosity])
        timesteps_sim = Jutul.report_timesteps(reports)
        cell_volume = prod(S.model.d)
        n_cells = length(ϕ)
        dϕ = zeros(T, n_cells)
        epsilon = T(1e-4)
        
        # Compute base objective
        G_base = zero(T)
        for (step_ix, (s, dt)) in enumerate(zip(states, timesteps_sim))
            step_info = Dict{Symbol, Any}(:step => step_ix, :ministep => 1)
            G_base += weighted_objective(model, s, dt, step_info, forces)
        end
        
        # Compute gradient for each cell using central difference
        for i in 1:n_cells
            if ϕ[i] >= 1.0
                dϕ[i] = zero(T)
                continue
            end
            
            # Forward perturbation
            ϕ_plus = copy(ϕ)
            ϕ_plus[i] += epsilon
            
            parameters_plus = deepcopy(parameters)
            parameters_plus[:FluidVolume] = cell_volume .* ϕ_plus
            
            model_plus = deepcopy(model)
            model_plus.data_domain[:porosity] = ϕ_plus
            
            states_plus, _ = simulate(dict(state0_), model_plus, tstep, parameters = parameters_plus, forces = forces, info_level=-1, max_timestep_cuts = 1000)
            
            G_plus = zero(T)
            for (step_ix, (s, dt)) in enumerate(zip(states_plus, timesteps_sim))
                step_info = Dict{Symbol, Any}(:step => step_ix, :ministep => 1)
                G_plus += weighted_objective(model_plus, s, dt, step_info, forces)
            end
            
            # Backward perturbation
            ϕ_minus = copy(ϕ)
            ϕ_minus[i] -= epsilon
            
            parameters_minus = deepcopy(parameters)
            parameters_minus[:FluidVolume] = cell_volume .* ϕ_minus
            
            model_minus = deepcopy(model)
            model_minus.data_domain[:porosity] = ϕ_minus
            
            states_minus, _ = simulate(dict(state0_), model_minus, tstep, parameters = parameters_minus, forces = forces, info_level=-1, max_timestep_cuts = 1000)
            
            G_minus = zero(T)
            for (step_ix, (s, dt)) in enumerate(zip(states_minus, timesteps_sim))
                step_info = Dict{Symbol, Any}(:step => step_ix, :ministep => 1)
                G_minus += weighted_objective(model_minus, s, dt, step_info, forces)
            end
            
            # Central difference
            dϕ[i] = (G_plus - G_minus) / (T(2) * epsilon)
        end
        
        if info_level >= 0
            println("DEBUG: dLogT[1:5] = ", dLogTransmissibilities[1:min(5, length(dLogTransmissibilities))])
            println("DEBUG: dϕ[1:5] = ", dϕ[1:min(5, length(dϕ))])
        end
        
        return NoTangent(), dLogTransmissibilities, dϕ, NoTangent()
    end
    return output, pullback
end

"""
    loss_per_step(m, state, dt, step_no, forces, states_ref, reservoir_states=nothing)

Compute loss for a single time step by comparing state with reference state.
"""
function loss_per_step(m, state, dt, step_no, forces, states_ref, reservoir_states=nothing)
    step_ix = _safe_step_index(step_no, length(states_ref))
    if reservoir_states !== nothing
        @assert length(states_ref) == length(reservoir_states) "states_ref length $(length(states_ref)) != reservoir_states length $(length(reservoir_states))"
    end
    state_ref = states_ref[step_ix]
    val = state[:Reservoir][:Saturations]
    val2 = state[:Reservoir][:Pressure]
    ref = state_ref[:Reservoir][:Saturations]
    ref2 = state_ref[:Reservoir][:Pressure]
    return inner_mismatch(val, ref, val2, ref2)
end

"""
    loss_per_step_simple(m, state, dt, step_no, forces, states_ref)

Compute loss for a single time step (simple model version).
"""
function loss_per_step_simple(m, state, dt, step_no, forces, states_ref)
    step_ix = _safe_step_index(step_no, length(states_ref))
    state_ref = states_ref[step_ix]
    val = state[:Saturations]
    val2 = state[:Pressure]
    ref = state_ref[:Saturations]
    ref2 = state_ref[:Pressure]
    return inner_mismatch(val, ref, val2, ref2)
end

"""
    inner_mismatch(val, ref, val2, ref2)

Compute squared mismatch between values and references for both saturation and pressure.
"""
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
    # Initialize fields required by newer Jutul versions
    data[:include_state0] = false
    data[:best_obj] = Inf
    data[:best_x] = copy(x0)  # Initialize with x0

    sim = Simulator(case)
    if isnothing(config)
        config = simulator_config(sim; info_level = -1, kwarg...)
    elseif !verbose
        config[:info_level] = -1
        config[:end_report] = false
    end
    data[:sim] = sim
    data[:sim_config] = config
    # Set include_state0 in data dict for objective_opt! (required in newer Jutul versions)
    if !haskey(data, :include_state0)
        data[:include_state0] = false
    end

    if grad_type == :adjoint
        # n_objective = nothing means single scalar objective (returns vector)
        # n_objective = 1 means single objective but returns matrix format
        # For compatibility, use nothing to get vector format
        adj_storage = setup_adjoint_storage(model, state0 = state0,
                                                   parameters = parameters,
                                                   targets = targets,
                                                   use_sparsity = use_sparsity,
                                                   param_obj = param_obj,
                                                   n_objective = nothing)
        data[:adjoint_storage] = adj_storage
        # grad_adj will be filled by solve_adjoint_sensitivities! in gradient_opt!
        # Initialize to match storage[:dparam] format (which uses gradient_vec_or_mat)
        # When n_objective = 1, gradient_vec_or_mat returns a column vector (matrix format)
        # We need to match this format, but transfer_gradient! expects a vector
        # So we initialize as vector and let solve_adjoint_sensitivities! handle the format
        grad_adj = zeros(adj_storage.n)
    else
        grad_adj = similar(x0)
    end
    data[:case] = case
    data[:grad_adj] = grad_adj
    data[:mapper] = mapper
    data[:G] = G
    data[:targets] = targets
    data[:config] = opt_cfg
    data[:last_obj] = Inf
    data[:x_hash] = hash(x0)
    data[:states] = precomputed_states
    data[:reports] = reports
    data[:dt_current] = case.dt
    F = x -> objective_opt!(x, data, print)
    dF = (dFdx, x) -> gradient_opt!(dFdx, x, data)
    F_and_dF = (F, dFdx, x) -> objective_and_gradient_opt!(F, dFdx, x, data, print)
    return (F! = F, dF! = dF, F_and_dF! = F_and_dF, x0 = x0, limits = lims, data = data)
end

