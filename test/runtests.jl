using MTk2JuMP
using JuMP
using Ipopt
using InteractiveUtils
using Test

const CART_POLE_DIR = joinpath(@__DIR__, "..", "examples", "cart_pole")
const HEV_DIR = joinpath(@__DIR__, "..", "examples", "HEV")
const ROCKET_LANDING_DIR = joinpath(@__DIR__, "..", "examples", "rocket_landing")
const IFB = MTk2JuMP.IF.build
const IFLG = MTk2JuMP.IF.LGPoly
const MOI = JuMP.MOI
include(joinpath(CART_POLE_DIR, "model.jl"))
include(joinpath(HEV_DIR, "model.jl"))
include(joinpath(ROCKET_LANDING_DIR, "model.jl"))
const CART_POLE_MODEL = CartPoleModel
const HEV_MODEL = HEVModel
const ROCKET_LANDING_MODEL = RocketLandingModel
const CART_POLE_CONFIG = include(joinpath(CART_POLE_DIR, "Config.jl"))
const CART_POLE_SB = include(joinpath(CART_POLE_DIR, "ScalesBounds.jl"))
const HEV_CONFIG = include(joinpath(HEV_DIR, "Config.jl"))
const HEV_SB = include(joinpath(HEV_DIR, "ScalesBounds.jl"))
const ROCKET_LANDING_CONFIG = include(joinpath(ROCKET_LANDING_DIR, "Config.jl"))
const ROCKET_LANDING_SB = include(joinpath(ROCKET_LANDING_DIR, "ScalesBounds.jl"))

function prepare_cart_pole_pipeline(n::Int; int_method::Symbol=:Coll, coll_poly::Symbol=:LGR, coll_order::Int=3)
    ocpi = IFB.OCPInterface()
    ocpi.model = Model(Ipopt.Optimizer; add_bridges=false)

    config = deepcopy(CART_POLE_CONFIG)
    config.Discretization.N = n
    config.Integration.int_method = int_method
    config.Integration.Collocation.poly = coll_poly
    config.Integration.Collocation.order = coll_order

    IFB.setup_problem!(ocpi, config)
    ocpi.sys = CART_POLE_MODEL.cart_pole()
    IFB.get_ode_architecture!(ocpi; verbose=false)
    IFB.set_bounds!(ocpi, CART_POLE_SB.Bounds)
    IFB.set_scales!(ocpi, CART_POLE_SB.Scales)

    x0 = zeros(Float64, ocpi.meta.nx, ocpi.settings.Discretization.N)
    u0 = zeros(Float64, ocpi.meta.nu, ocpi.settings.Discretization.N - 1)
    IFB.set_opt_vars!(ocpi; x0=x0, u0=u0)
    IFB.set_y_expr!(ocpi)
    IFB.set_gDyn!(ocpi; f_scale=ocpi.settings.di)
    IFB.set_sparse_param_closure!(ocpi)
    IFB.set_control_der!(ocpi; dt=[ocpi.settings.di])
    IFB.set_control_acc!(ocpi; dt=[ocpi.settings.di])

    x_idx = findfirst(==("x(t)"), ocpi.meta.x_names)
    v_idx = findfirst(==("v(t)"), ocpi.meta.x_names)
    θ_idx = findfirst(==("θ(t)"), ocpi.meta.x_names)
    ω_idx = findfirst(==("ω(t)"), ocpi.meta.x_names)
    a_idx = findfirst(==("a(t)"), ocpi.meta.x_names)

    @constraint(ocpi.model, ocpi.vars.x[x_idx, 1] == 0.0)
    @constraint(ocpi.model, ocpi.vars.x[x_idx, end] == 0.0)
    @constraint(ocpi.model, ocpi.vars.x[v_idx, 1] == 0.0)
    @constraint(ocpi.model, ocpi.vars.x[v_idx, end] == 0.0)
    @constraint(ocpi.model, ocpi.vars.x[θ_idx, 1] == 0.0)
    @constraint(ocpi.model, ocpi.vars.x[θ_idx, end] == pi)
    @constraint(ocpi.model, ocpi.vars.x[ω_idx, 1] == 0.0)
    @constraint(ocpi.model, ocpi.vars.x[ω_idx, end] == 0.0)
    @constraint(ocpi.model, ocpi.vars.x[a_idx, 1] == 0.0)
    @constraint(ocpi.model, ocpi.vars.x[a_idx, end] == 0.0)

    P_expr = IFB.get_yi_by_name(ocpi, :P)
    @objective(ocpi.model, Min, sum(P_expr[i]^2 for i in 1:(ocpi.settings.Discretization.N - 1)))

    return ocpi
end

function build_cart_pole_pipeline(n::Int)
    ocpi = prepare_cart_pole_pipeline(n)
    optimize!(ocpi.model)
    IFB.collect_all_results!(ocpi)
    IFB.collect_solver_info!(ocpi)

    return ocpi
end

function prepare_rocket_landing_pipeline(n::Int; int_method::Symbol=:Coll, coll_poly::Symbol=:LGR, coll_order::Int=3)
    ocpi = IFB.OCPInterface()
    ocpi.model = Model(Ipopt.Optimizer; add_bridges=false)

    config = deepcopy(ROCKET_LANDING_CONFIG)
    config.Discretization.N = n
    config.Integration.int_method = int_method
    config.Integration.Collocation.poly = coll_poly
    config.Integration.Collocation.order = coll_order

    IFB.setup_problem!(ocpi, config)
    ocpi.sys = ROCKET_LANDING_MODEL.rocket_landing()
    IFB.get_ode_architecture!(ocpi; verbose=false)
    IFB.set_bounds!(ocpi, ROCKET_LANDING_SB.Bounds)
    IFB.set_scales!(ocpi, ROCKET_LANDING_SB.Scales)

    x0 = zeros(Float64, ocpi.meta.nx, ocpi.settings.Discretization.N)
    u0 = zeros(Float64, ocpi.meta.nu, ocpi.settings.Discretization.N - 1)
    IFB.set_opt_vars!(ocpi; x0=x0, u0=u0)
    IFB.set_y_expr!(ocpi)
    IFB.set_gDyn!(ocpi; f_scale=ocpi.settings.di)
    IFB.set_sparse_param_closure!(ocpi)
    IFB.set_control_der!(ocpi; dt=[ocpi.settings.di])
    IFB.set_control_acc!(ocpi; dt=[ocpi.settings.di])

    x_idx = findfirst(==("x(t)"), ocpi.meta.x_names)
    y_idx = findfirst(==("y(t)"), ocpi.meta.x_names)
    vx_idx = findfirst(==("v_x(t)"), ocpi.meta.x_names)
    vy_idx = findfirst(==("v_y(t)"), ocpi.meta.x_names)
    θ_idx = findfirst(==("θ(t)"), ocpi.meta.x_names)
    ω_idx = findfirst(==("ω(t)"), ocpi.meta.x_names)

    @constraint(ocpi.model, ocpi.vars.x[x_idx, 1] == 100.0)
    @constraint(ocpi.model, ocpi.vars.x[y_idx, 1] == 100.0)
    @constraint(ocpi.model, ocpi.vars.x[vx_idx, 1] == 0.0)
    @constraint(ocpi.model, ocpi.vars.x[vx_idx, end] == 0.0)
    @constraint(ocpi.model, ocpi.vars.x[vy_idx, 1] == -5.0)
    @constraint(ocpi.model, ocpi.vars.x[vy_idx, end] == 0.0)
    @constraint(ocpi.model, ocpi.vars.x[θ_idx, 1] == pi / 6)
    @constraint(ocpi.model, ocpi.vars.x[θ_idx, end] == 0.0)
    @constraint(ocpi.model, ocpi.vars.x[ω_idx, 1] == 0.2)
    @constraint(ocpi.model, ocpi.vars.x[ω_idx, end] == 0.0)

    y_base_expr = IFB.get_yi_by_name(ocpi, :y_base)
    x_base_expr = IFB.get_yi_by_name(ocpi, :x_base)
    @constraint(ocpi.model, [i in 1:(ocpi.settings.N - 1)], y_base_expr[i] >= 0.0)
    @constraint(ocpi.model, y_base_expr[end] == 0.0)
    @constraint(ocpi.model, x_base_expr[end] == 0.0)

    obj = sum(ocpi.vars_n.u[i, j]^2 for i in 1:ocpi.meta.nu, j in 1:(ocpi.settings.N - 1))
    reg = sum((ocpi.vars_n.x_col[θ_idx, j, k])^2 for j in 1:(ocpi.settings.N - 1), k in 1:ocpi.settings.Coll_set.order)
    @objective(ocpi.model, Min, obj + 0.5 * reg)

    return ocpi
end

function prepare_hev_pipeline(n::Int; max_iter::Int=5)
    ocpi = IFB.OCPInterface()
    ocpi.model = Model(Ipopt.Optimizer; add_bridges=false)

    config = deepcopy(HEV_CONFIG)
    config.Discretization.N = n
    config.Ipopt.max_iter = max_iter

    IFB.setup_problem!(ocpi, config)
    IFB.load_LUT2d!(ocpi, joinpath(HEV_DIR, "bsfc_lut.json"), :bsfc_lut; x1_name="omega", x2_name="torque", y_name="bsfc")
    IFB.load_LUT1d!(ocpi, joinpath(HEV_DIR, "r0_lut.json"), :r0_lut; x_name="SOC", y_name="r0")
    IFB.load_LUT1d!(ocpi, joinpath(HEV_DIR, "voc_lut.json"), :voc_lut; x_name="SOC", y_name="voc")
    ocpi.sys = HEV_MODEL.hev(ocpi.LUTs1D, ocpi.LUTs2D)
    IFB.get_ode_architecture!(ocpi; verbose=false)
    IFB.set_bounds!(ocpi, HEV_SB.Bounds)
    IFB.set_scales!(ocpi, HEV_SB.Scales)

    x0 = zeros(Float64, ocpi.meta.nx, ocpi.settings.Discretization.N)
    u0 = zeros(Float64, ocpi.meta.nu, ocpi.settings.Discretization.N - 1)
    IFB.set_opt_vars!(ocpi; x0=x0, u0=u0)
    IFB.set_y_expr!(ocpi)
    IFB.set_gDyn!(ocpi; f_scale=ocpi.settings.di)
    IFB.set_sparse_param_closure!(ocpi)
    IFB.set_control_der!(ocpi; dt=[ocpi.settings.di])
    IFB.set_control_acc!(ocpi; dt=[ocpi.settings.di])

    x_idx = findfirst(==("x(t)"), ocpi.meta.x_names)
    SOC_idx = findfirst(==("SOC(t)"), ocpi.meta.x_names)
    m_fuel_idx = findfirst(==("m_fuel(t)"), ocpi.meta.x_names)
    v_idx = findfirst(==("v(t)"), ocpi.meta.x_names)

    @constraint(ocpi.model, ocpi.vars.x[x_idx, 1] == 0.0)
    @constraint(ocpi.model, ocpi.vars.x[x_idx, end] == 100.0)
    @constraint(ocpi.model, ocpi.vars.x[SOC_idx, 1] == 1.0)
    @constraint(ocpi.model, ocpi.vars.x[m_fuel_idx, 1] == 0.0)
    @constraint(ocpi.model, ocpi.vars.x[v_idx, 1] == 2.0)

    reg = sum(ocpi.vars_n.u[i, j]^2 for i in 1:ocpi.meta.nu, j in 1:(ocpi.settings.N - 1))
    @objective(ocpi.model, Min, 0.01 * reg)

    return ocpi
end

@testset "MTk2JuMP.jl" begin
    @testset "problem config application" begin
        config = IFB.OCPSettings(
            Ipopt = IFB.IpoptConfig(
                tol = 1e-6,
                max_iter = 10_000,
                hsllib_path = "",
                linear_solver = "ma27",
                warm_start = true,
            ),
            Discretization = IFB.DiscretizationConfig(
                N = 12,
                tspan = (0.0, 6.0),
            ),
            Integration = IFB.IntegrationConfig(
                int_method = :Coll,
                Collocation = IFB.CollocationConfig(
                    poly = :LGR,
                    order = 3,
                ),
            ),
            Params = IFB.ParamsConfig(mode = :sparse),
        )

        ocpi = IFB.OCPInterface()
    IFB.setup_problem!(ocpi, config)

        @test ocpi.settings.Discretization.N == 12
        @test ocpi.settings.Discretization.tspan == (0.0, 6.0)
        @test isapprox(ocpi.settings.di, 6.0 / 11.0)
        @test ocpi.settings.Integration.int_method == :Coll
        @test ocpi.settings.Integration.Collocation.order == 3
        @test ocpi.settings.Integration.Collocation.poly == :LGR
        @test ocpi.settings.Params.mode == :sparse
    end

    @testset "invalid problem config" begin
        bad_config = IFB.OCPSettings(
            Ipopt = IFB.IpoptConfig(
                tol = 1e-6,
                max_iter = 10_000,
                hsllib_path = "",
                linear_solver = "ma27",
                warm_start = true,
            ),
            Discretization = IFB.DiscretizationConfig(
                N = 1,
                tspan = (0.0, 1.0),
            ),
            Integration = IFB.IntegrationConfig(
                int_method = :Coll,
                Collocation = IFB.CollocationConfig(),
            ),
            Params = IFB.ParamsConfig(mode = :single),
        )

        ocpi = IFB.OCPInterface()
        @test_throws ErrorException IFB.validate_settings!(bad_config)
    end

    @testset "invalid derived settings" begin
        config = IFB.OCPSettings(
            Ipopt = IFB.IpoptConfig(
                tol = 1e-6,
                max_iter = 10_000,
                hsllib_path = "",
                linear_solver = "ma27",
                warm_start = true,
            ),
            Discretization = IFB.DiscretizationConfig(
                N = 12,
                tspan = (0.0, 6.0),
            ),
            Integration = IFB.IntegrationConfig(
                int_method = :Coll,
                Collocation = IFB.CollocationConfig(
                    poly = :LGR,
                    order = 0,
                ),
            ),
            Params = IFB.ParamsConfig(mode = :sparse),
        )

        ocpi = IFB.OCPInterface()
        @test_throws ErrorException IFB.validate_settings!(config)
        @test_throws ErrorException IFB.setup_problem!(ocpi, config)
    end

    @testset "solver info collection" begin
        ocpi = IFB.OCPInterface()
        ocpi.model = Model(Ipopt.Optimizer; add_bridges=false)

        @variable(ocpi.model, x >= 0)
        @objective(ocpi.model, Min, (x - 2)^2)
        optimize!(ocpi.model)

        info = IFB.collect_solver_info!(ocpi)

        @test info.termination_status == JuMP.termination_status(ocpi.model)
        @test isfinite(info.solver_time)
        @test info.iterations ≥ 0
        @test isapprox(info.objective_value, 0.0; atol=1e-6)
        @test isfinite(info.dual_value)
    end

    @testset "collocation method coverage" begin
        for method in (:LG, :LGL)
            nodes, C, B, D = IFLG.get_collocation_system(4; method=method)

            @test length(nodes) == 5
            @test size(C) == (5, 5)
            @test length(B) == 5
            @test length(D) == 5
            @test isapprox(nodes[1], 0.0; atol=1e-12)
            @test all(isfinite, nodes)
            @test all(isfinite, C)
            @test all(isfinite, B)
            @test all(isfinite, D)
        end

        lgl_nodes, _, _, _ = IFLG.get_collocation_system(4; method=:LGL)
        @test isapprox(lgl_nodes[end], 1.0; atol=1e-12)

        lg_nodes, _, _, _ = IFLG.get_collocation_system(4; method=:LG)
        @test lg_nodes[end] < 1.0
    end

    @testset "EE integration method coverage" begin
        ocpi = prepare_cart_pole_pipeline(80; int_method=:EE)
        optimize!(ocpi.model)

        @test JuMP.termination_status(ocpi.model) in (MOI.LOCALLY_SOLVED, MOI.OPTIMAL)
        @test size(ocpi.vars.x) == (ocpi.meta.nx, ocpi.settings.N)
        @test size(ocpi.vars.u) == (ocpi.meta.nu, ocpi.settings.N - 1)
        @test isempty(ocpi.vars.x_col)
    end

    @testset "LG and LGL integration coverage" begin
        for method in (:LG, :LGL)
            ocpi = prepare_cart_pole_pipeline(70; int_method=:Coll, coll_poly=method)
            optimize!(ocpi.model)

            @test JuMP.termination_status(ocpi.model) in (MOI.LOCALLY_SOLVED, MOI.OPTIMAL)
            @test ocpi.settings.Coll_set.poly == method
            @test !isempty(ocpi.vars.x_col)

            IFB.collect_all_results!(ocpi)
            @test haskey(ocpi.res.x, :x)
            @test haskey(ocpi.res.y, :P)
        end
    end

    @testset "rocket landing example integration" begin
        ocpi = prepare_rocket_landing_pipeline(45; int_method=:Coll, coll_poly=:LGR)
        optimize!(ocpi.model)
        IFB.collect_all_results!(ocpi)
        IFB.collect_solver_info!(ocpi)

        @test JuMP.termination_status(ocpi.model) in (MOI.LOCALLY_SOLVED, MOI.OPTIMAL)
        @test haskey(ocpi.res.x, :x)
        @test haskey(ocpi.res.x, :y)
        @test haskey(ocpi.res.u, :F_m)
        @test haskey(ocpi.res.y, :x_base)
        @test haskey(ocpi.res.y, :y_base)
        @test isfinite(ocpi.solver_info.objective_value)
    end

    @testset "HEV smooth interpolation integration" begin
        ocpi = prepare_hev_pipeline(30; max_iter=5)

        voc_lut = ocpi.LUTs1D[:voc_lut]
        r0_lut = ocpi.LUTs1D[:r0_lut]
        bsfc_lut = ocpi.LUTs2D[:bsfc_lut]

        voc_mid = cld(length(voc_lut.x), 2)
        r0_mid = cld(length(r0_lut.x), 2)
        x1_mid = cld(length(bsfc_lut.x1), 2)
        x2_mid = cld(length(bsfc_lut.x2), 2)

        @test isapprox(voc_lut.itp(voc_lut.x[voc_mid]), voc_lut.y[voc_mid]; atol=1e-6, rtol=1e-6)
        @test isapprox(r0_lut.itp(r0_lut.x[r0_mid]), r0_lut.y[r0_mid]; atol=1e-6, rtol=1e-6)
        @test isapprox(bsfc_lut.itp(bsfc_lut.x1[x1_mid], bsfc_lut.x2[x2_mid]), bsfc_lut.y[x1_mid, x2_mid]; atol=1e-5, rtol=1e-5)

        optimize!(ocpi.model)

        @test JuMP.termination_status(ocpi.model) in (MOI.LOCALLY_SOLVED, MOI.OPTIMAL, MOI.ITERATION_LIMIT)
        @test isfinite(JuMP.objective_value(ocpi.model))
    end

    @testset "SmoothInterpolations accuracy against ground truth" begin
        @testset "1D smooth interpolator accuracy" begin
            f1d(x) = sin(1.2 * x) + 0.1 * x^2

            x_train = collect(range(-1.5, 1.5; length=15))
            y_train = f1d.(x_train)
            itp1d = MTk2JuMP.IF.SmoothInterpolations.build_smooth_rbf1d(x_train, y_train)

            train_err = maximum(abs.([itp1d(xi) - yi for (xi, yi) in zip(x_train, y_train)]))
            @test train_err < 1e-9

            x_test = collect(range(-1.4, 1.4; length=29))
            y_test = f1d.(x_test)
            test_err = maximum(abs.([itp1d(xi) - yi for (xi, yi) in zip(x_test, y_test)]))
            @test test_err < 0.09
        end

        @testset "2D smooth interpolator accuracy" begin
            f2d(x, y) = sin(1.1 * x) * cos(0.9 * y) + 0.05 * x^2 - 0.03 * y

            x1_train = collect(range(-1.0, 1.0; length=9))
            x2_train = collect(range(-0.8, 0.8; length=8))
            z_train = [f2d(x1, x2) for x1 in x1_train, x2 in x2_train]
            itp2d = MTk2JuMP.IF.SmoothInterpolations.build_smooth_rbf2d(x1_train, x2_train, z_train)

            train_err = maximum(abs.([
                itp2d(x1, x2) - z_train[i, j] for (i, x1) in enumerate(x1_train), (j, x2) in enumerate(x2_train)
            ]))
            @test train_err < 1e-9

            x1_test = collect(range(-0.9, 0.9; length=7))
            x2_test = collect(range(-0.7, 0.7; length=6))
            test_err = maximum(abs.([
                itp2d(x1, x2) - f2d(x1, x2) for x1 in x1_test, x2 in x2_test
            ]))
            @test test_err < 0.11
        end
    end

    @testset "cart pole pipeline" begin
        @code_warntype prepare_cart_pole_pipeline(200)
        ocpi = build_cart_pole_pipeline(200)
        @test JuMP.termination_status(ocpi.model) in (MOI.LOCALLY_SOLVED, MOI.OPTIMAL)
        @test haskey(ocpi.res.x, :x)
        @test haskey(ocpi.res.u, :F)
        @test haskey(ocpi.res.y, :P)
        @test haskey(ocpi.duals.gDyn, :x)
        @test isfinite(ocpi.solver_info.dual_value)
    end
end
