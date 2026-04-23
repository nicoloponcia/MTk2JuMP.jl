function validate_settings!(NLPConfig::OCPSettings_)
    NLPConfig.Discretization.N ≥ 2 || error("Discretization.N must be at least 2.")
    isfinite(NLPConfig.Discretization.tspan[1]) || error("Discretization.tspan start must be finite.")
    isfinite(NLPConfig.Discretization.tspan[2]) || error("Discretization.tspan stop must be finite.")
    NLPConfig.Discretization.tspan[2] > NLPConfig.Discretization.tspan[1] || error("Discretization.tspan must increase from start to stop.")
    NLPConfig.Integration.int_method ∈ (:Coll, :EE) || error("Unsupported integration method: $(NLPConfig.Integration.int_method). Supported methods are :Coll and :EE.")
    NLPConfig.Params.mode ∈ (:single, :sparse) || error("Unsupported parameter mode: $(NLPConfig.Params.mode). Supported modes are :single and :sparse.")

    if NLPConfig.Integration.int_method == :Coll
        NLPConfig.Integration.Collocation.order > 0 || error("Collocation.order must be positive.")
        NLPConfig.Integration.Collocation.poly ∈ (:LG, :LGR, :LGL) || error("Unsupported collocation polynomial: $(NLPConfig.Integration.Collocation.poly). Supported methods are :LG, :LGR, and :LGL.")
    end

    return NLPConfig
end

function _validate_optimizer_config!(NLPConfig)
    NLPConfig.Ipopt.max_iter > 0 || error("Ipopt.max_iter must be positive.")
    NLPConfig.Ipopt.tol > 0 || error("Ipopt.tol must be positive.")
end

function set_settings!(OCPI::OCPInterface_, NLPConfig::OCPSettings_)
    OCPI.settings = validate_settings!(NLPConfig)
    if OCPI.settings.Integration.int_method == :Coll
        LGPoly.set_coll_method!(OCPI.settings.Integration.Collocation)
    end
    return OCPI
end

function setup_optimizer!(model::Model, NLPConfig::OCPSettings_)
    _validate_optimizer_config!(NLPConfig)

    set_optimizer_attribute(model, "tol", NLPConfig.Ipopt.tol)
    set_attribute(model, "hsllib", NLPConfig.Ipopt.hsllib_path)
    set_attribute(model, "linear_solver", NLPConfig.Ipopt.linear_solver)
    set_attribute(model, "max_iter", NLPConfig.Ipopt.max_iter)
    set_string_names_on_creation(model, false)
    set_attribute(model, "warm_start_init_point", NLPConfig.Ipopt.warm_start ? "yes" : "no")
    
    set_optimizer_attribute(model, "print_timing_statistics", "yes")
    set_optimizer_attribute(model, "print_level", 5)
    
    set_attribute(model, "hessian_approximation", "exact")
    # set_optimizer_attribute(model, "hessian_approximation", "limited-memory")
    return model
end



function setup_problem!(OCPI::OCPInterface_, NLPConfig::OCPSettings_)
    OCPI.settings = validate_settings!(NLPConfig)
    validate_settings!(OCPI.settings)
    if OCPI.settings.Integration.int_method == :Coll
        LGPoly.set_coll_method!(OCPI.settings.Integration.Collocation)
    end
    setup_optimizer!(OCPI.model, OCPI.settings)
    return OCPI
end




