function get_ode_architecture!(OCPI::OCPInterface_; verbose::Bool=false)
    t = ModelingToolkit.get_iv(OCPI.sys)

    vars   = ModelingToolkit.unknowns(OCPI.sys)
    params = ModelingToolkit.parameters(OCPI.sys)
    obs    = ModelingToolkit.observables(OCPI.sys)
    defs   = ModelingToolkit.initial_conditions(OCPI.sys)


    inputs = filter(params) do p
        istree(p) && any(isequal(t), arguments(p))
    end

    opt_params = filter(params) do p
        is_static = !(istree(p) && any(isequal(t), arguments(p)))
        has_no_default = !haskey(defs, p)
        return is_static && has_no_default
    end

    OCPI.meta.nx = length(vars)
    OCPI.meta.nu = length(inputs)
    OCPI.meta.np = length(opt_params)
    OCPI.meta.ny = length(obs)

    OCPI.meta.x_names = [string(v) for v in vars]
    OCPI.meta.u_names = [string(u) for u in inputs]
    OCPI.meta.p_names = [string(p) for p in opt_params]
    OCPI.meta.y_names = [string(y) for y in obs]

    OCPI.decision_vars = vcat(vars, inputs, opt_params)

    if verbose
        println("OCP Architecture Summary:")
        println("  Number of state variables (nx): ", OCPI.meta.nx)
        println("  Number of control inputs (nu): ", OCPI.meta.nu)
        println("  Number of optimization parameters (np): ", OCPI.meta.np)
        println("  State variable names: ", OCPI.meta.x_names)
        println("  Control input names: ", OCPI.meta.u_names)
        println("  Optimization parameter names: ", OCPI.meta.p_names)
    end

end
