using MTk2JuMP
using JuMP
using Ipopt
using Test

@testset "MTk2JuMP.jl" begin
    IFB = MTk2JuMP.IF.build

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
        @test haskey(info.primals, :decision_vars)
        @test length(info.primals[:decision_vars]) == 1
        @test isapprox(info.primals[:decision_vars][1], 2.0; atol=1e-6)
    end
end
