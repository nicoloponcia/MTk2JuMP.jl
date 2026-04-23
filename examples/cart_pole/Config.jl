function get_nlp_config()
    Ipopt = (
        tol = 1e-6,
        max_iter = 10000,
        hsllib_path = "C:\\ProgramData\\HSL\\bin\\libhsl.dll",
        linear_solver = "ma27",
        warm_start = true,
    )

    Discretization = (
        N = 50,
        tspan = (0.0, 5.0),
    )

    Collocation = (
        poly = :LGR,
        order = 3,
    )


    Integration = (
        int_method = :Coll,
        Collocation = Collocation,
    )

    Params = (
        mode = :sparse, # :sparse or :single
    )

    # Derived scales
    NLPConfig = (
        Ipopt = Ipopt,
        Discretization = Discretization,
        Integration = Integration,
        Params = Params
    )
return NLPConfig
end

NLPConfig = get_nlp_config()