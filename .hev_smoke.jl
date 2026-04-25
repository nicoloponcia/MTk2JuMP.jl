using MTk2JuMP
using JuMP
using Ipopt

const HEV_DIR = joinpath(@__DIR__, "examples", "HEV")
include(joinpath(HEV_DIR, "model.jl"))
config = deepcopy(include(joinpath(HEV_DIR, "Config.jl")))
sb = include(joinpath(HEV_DIR, "ScalesBounds.jl"))
config.Discretization.N = 8
config.Ipopt.max_iter = 1

ocpi = MTk2JuMP.IF.build.OCPInterface()
ocpi.model = Model(Ipopt.Optimizer; add_bridges=false)
MTk2JuMP.IF.build.setup_problem!(ocpi, config)
MTk2JuMP.IF.build.load_LUT2d!(ocpi, joinpath(HEV_DIR, "bsfc_lut.json"), :bsfc_lut; x1_name="omega", x2_name="torque", y_name="bsfc")
MTk2JuMP.IF.build.load_LUT1d!(ocpi, joinpath(HEV_DIR, "r0_lut.json"), :r0_lut; x_name="SOC", y_name="r0")
MTk2JuMP.IF.build.load_LUT1d!(ocpi, joinpath(HEV_DIR, "voc_lut.json"), :voc_lut; x_name="SOC", y_name="voc")
ocpi.sys = HEVModel.hev(ocpi.LUTs1D, ocpi.LUTs2D)
MTk2JuMP.IF.build.get_ode_architecture!(ocpi; verbose=false)
MTk2JuMP.IF.build.set_bounds!(ocpi, sb.Bounds)
MTk2JuMP.IF.build.set_scales!(ocpi, sb.Scales)
x0 = zeros(Float64, ocpi.meta.nx, ocpi.settings.Discretization.N)
u0 = zeros(Float64, ocpi.meta.nu, ocpi.settings.Discretization.N - 1)
MTk2JuMP.IF.build.set_opt_vars!(ocpi; x0=x0, u0=u0)
MTk2JuMP.IF.build.set_y_expr!(ocpi)
MTk2JuMP.IF.build.set_gDyn!(ocpi; f_scale=ocpi.settings.di)
MTk2JuMP.IF.build.set_sparse_param_closure!(ocpi)
MTk2JuMP.IF.build.set_control_der!(ocpi; dt=[ocpi.settings.di])
MTk2JuMP.IF.build.set_control_acc!(ocpi; dt=[ocpi.settings.di])

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

optimize!(ocpi.model)
println(termination_status(ocpi.model))
println(isfinite(objective_value(ocpi.model)))
