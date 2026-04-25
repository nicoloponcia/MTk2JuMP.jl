function get_nlp_config()
    return MTk2JuMP.IF.build.OCPSettings(
        Ipopt = MTk2JuMP.IF.build.IpoptConfig(
        tol = 1e-6,
        max_iter = 10000,
        # hsllib_path = "C:\\ProgramData\\HSL\\bin\\libhsl.dll",
        hsllib_path = "",
        # linear_solver = "ma27",
        linear_solver = "mumps",
        warm_start = true,
        ),
        Discretization = MTk2JuMP.IF.build.DiscretizationConfig(
        N = 100,
        tspan = (0.0, 5.0),
        ),
        Integration = MTk2JuMP.IF.build.IntegrationConfig(
        int_method = :Coll,
        Collocation = MTk2JuMP.IF.build.CollocationConfig(
        poly = :LGR,
        order = 3,
        ),
        ),
        Params = MTk2JuMP.IF.build.ParamsConfig(
        mode = :sparse, # :sparse or :single
        )
    )
end

NLPConfig = get_nlp_config()