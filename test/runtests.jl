using MTk2JuMP
using JuMP
using Ipopt
using InteractiveUtils
using Test

const CART_POLE_DIR = joinpath(@__DIR__, "..", "examples", "cart_pole")
const IFB = MTk2JuMP.IF.build
const MOI = JuMP.MOI
include(joinpath(CART_POLE_DIR, "model.jl"))
const CART_POLE_MODEL = CartPoleModel
const CART_POLE_CONFIG = include(joinpath(CART_POLE_DIR, "Config.jl"))
const CART_POLE_SB = include(joinpath(CART_POLE_DIR, "ScalesBounds.jl"))

function prepare_cart_pole_pipeline(n::Int)
    ocpi = IFB.OCPInterface()
    ocpi.model = Model(Ipopt.Optimizer; add_bridges=false)

    config = deepcopy(CART_POLE_CONFIG)
    config.Discretization.N = n

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

    @testset "cart pole pipeline" begin
        @inferred prepare_cart_pole_pipeline(200)
        ocpi = build_cart_pole_pipeline(200)
        @test JuMP.termination_status(ocpi.model) in (MOI.LOCALLY_SOLVED, MOI.OPTIMAL)
        @test haskey(ocpi.res.x, :x)
        @test haskey(ocpi.res.u, :F)
        @test haskey(ocpi.res.y, :P)
        @test haskey(ocpi.duals.gDyn, :x)
        @test isfinite(ocpi.solver_info.dual_value)
    end
end
