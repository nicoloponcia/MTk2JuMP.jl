using MTk2JuMP
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
        ocpi.settings = config
        IFB.setup_problem!(ocpi, ocpi.model)

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
        ocpi.settings = config
        @test_throws ErrorException IFB.validate_settings!(config)
        @test_throws ErrorException IFB.setup_problem!(ocpi, ocpi.model)
    end
end
