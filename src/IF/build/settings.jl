function set_settings!(OCPI::OCPInterface_, NLPConfig::NamedTuple)
    OCPI.settings.N = NLPConfig.Discretization.N
    OCPI.settings.tspan = NLPConfig.Discretization.tspan
    OCPI.settings.di = (NLPConfig.Discretization.tspan[2] - NLPConfig.Discretization.tspan[1])/(NLPConfig.Discretization.N-1)

    OCPI.settings.int_method = NLPConfig.Integration.int_method

    if OCPI.settings.int_method == :Coll
        if !haskey(NLPConfig.Integration, :Collocation)
            error("Collocation settings must be provided for Coll integration method.")
        end
        OCPI.settings.Coll_set.order = NLPConfig.Integration.Collocation.order
        OCPI.settings.Coll_set.poly = NLPConfig.Integration.Collocation.poly
        LGPoly.set_coll_method!(OCPI.settings.Coll_set)
    elseif OCPI.settings.int_method == :EE
    else
        error("Unsupported integration method: $(OCPI.settings.int_method). Supported methods are :LGR and :EE.")
    end

    OCPI.settings.Params_mode = NLPConfig.Params.mode

end

function setup_optimizer!(model::Model, NLPConfig::Any)
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




