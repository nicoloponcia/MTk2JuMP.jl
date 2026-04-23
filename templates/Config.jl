function get_nlp_config()
    # This template keeps the configuration local to the example or simulation
    # that includes it. Copy it per run and adjust only the values that are
    # specific to that problem instance.

    # Solver settings for the nonlinear optimizer.
    ipopt = MTk2JuMP.IF.build.IpoptConfig(
        tol = 1e-6,
        max_iter = 10_000,
        hsllib_path = "",
        linear_solver = "mumps",
        warm_start = true,
    )

    # Time discretization for the optimal control grid.
    # tspan also works in case of state-based discretization
    discretization = MTk2JuMP.IF.build.DiscretizationConfig(
        N = 200,
        tspan = (0.0, 5.0),
    )

    # Integration and collocation method used to transcribe the dynamics.
    integration = MTk2JuMP.IF.build.IntegrationConfig(
        int_method = :Coll,
        Collocation = MTk2JuMP.IF.build.CollocationConfig(
            poly = :LGR,
            order = 3,
        ),
    )

    # Parameter handling mode for the runtime model.
    params = MTk2JuMP.IF.build.ParamsConfig(
        mode = :sparse, # use :single when parameters are shared across the horizon
    )

    return MTk2JuMP.IF.build.OCPSettings(
        Ipopt = ipopt,
        Discretization = discretization,
        Integration = integration,
        Params = params,
    )
end

NLPConfig = get_nlp_config()