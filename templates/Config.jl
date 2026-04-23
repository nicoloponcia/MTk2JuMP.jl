
Ipopt = (
    tol = 1e-6,
    max_iter = 10000,
    hsllib_path = "C:\\ProgramData\\HSL\\bin\\libhsl.dll",
    linear_solver = "ma27",
    warm_start = true,
)


Collocation = (
    poly = :LGR,
    order = 3,
)


Integration = (
    int_method = :Coll,
    Collocation = Collocation,
)


RaceOptions = (
    StandingStart = false,
)

Params = (
    mode = :sparse, # :sparse or :single
)

# Derived scales
NLPConfig = (
    Ipopt = Ipopt,
    Integration = Integration,
    RaceOptions = RaceOptions,
    Params = Params
)


return NLPConfig