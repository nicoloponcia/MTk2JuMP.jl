using MTk2JuMP
using JuMP
using Ipopt
using Plots

# create OCPI object
OCPI = MTk2JuMP.IF.build.OCPInterface()

# instantiate JuMP model
OCPI.model = Model(Ipopt.Optimizer; add_bridges=false)

# load configs
include("Config.jl")
MTk2JuMP.IF.build.setup_problem!(OCPI, NLPConfig)

# create time vector (simple for time-based discretization)
t = LinRange(OCPI.settings.tspan[1], OCPI.settings.tspan[2], OCPI.settings.N)

# load LUTs
MTk2JuMP.IF.build.load_LUT2d!(OCPI, "examples/HEV/bsfc_lut.json", :bsfc_lut, x1_name="omega", x2_name="torque", y_name="bsfc")
MTk2JuMP.IF.build.load_LUT1d!(OCPI, "examples/HEV/r0_lut.json", :r0_lut, x_name="SOC", y_name="r0")
MTk2JuMP.IF.build.load_LUT1d!(OCPI, "examples/HEV/voc_lut.json", :voc_lut, x_name="SOC", y_name="voc")

# Load MTk model
include("model.jl")
OCPI.sys = HEVModel.hev(OCPI.LUTs1D, OCPI.LUTs2D)

# extract sys structure
MTk2JuMP.IF.build.get_ode_architecture!(OCPI; verbose=false)

# get and set bounds
SB = include("ScalesBounds.jl")
MTk2JuMP.IF.build.set_bounds!(OCPI, SB.Bounds)
MTk2JuMP.IF.build.set_scales!(OCPI, SB.Scales)

# create decision variables and set initial guess
x0 = zeros(Float64, OCPI.meta.nx, OCPI.settings.N)
u0 = zeros(Float64, OCPI.meta.nu, OCPI.settings.N-1)

MTk2JuMP.IF.build.set_opt_vars!(OCPI; x0=x0, u0=u0)

# create symbolic expression for auxiliary variables
MTk2JuMP.IF.build.set_y_expr!(OCPI)

# create dynamics constraints
f_scale = OCPI.settings.di
MTk2JuMP.IF.build.set_gDyn!(OCPI; f_scale=f_scale)

# closure sparse parameters constraint if needed
MTk2JuMP.IF.build.set_sparse_param_closure!(OCPI)

# set control derivatives and accelerations if needed
MTk2JuMP.IF.build.set_control_der!(OCPI; dt=[OCPI.settings.di])
MTk2JuMP.IF.build.set_control_acc!(OCPI; dt=[OCPI.settings.di])

# set boundary constraints
# e.g. below, could develop another environment to handle more complex OCP
x_idx = findfirst(x -> x == "x(t)", OCPI.meta.x_names)
@constraint(OCPI.model, OCPI.vars.x[x_idx, 1] == 0.0)

SOC_idx = findfirst(x -> x == "SOC(t)", OCPI.meta.x_names)
@constraint(OCPI.model, OCPI.vars.x[SOC_idx, 1] == 0.9)

m_fuel_idx = findfirst(x -> x == "m_fuel(t)", OCPI.meta.x_names)
@constraint(OCPI.model, OCPI.vars.x[m_fuel_idx, 1] == 0.0)

# constraint velocity to follow wltp cycle
include("wltp.jl")
vel_itp = MTk2JuMP.IF.SmoothInterpolations.build_smooth_rbf1d(wltp.time, wltp.velocity)

v_idx = findfirst(x -> x == "v(t)", OCPI.meta.x_names)
@constraint(OCPI.model, [j in 1:OCPI.settings.N], OCPI.vars.x[v_idx, j] == vel_itp(t[j]))


# define objective to minimize energy consumption
# P_expr = MTk2JuMP.IF.build.get_yi_by_name(OCPI, :P)
m_fuel_idx = findfirst(x -> x == "m_fuel(t)", OCPI.meta.x_names)
obj = OCPI.vars.x[m_fuel_idx, end]

# add regularization on control imputs
# reg = sum(OCPI.vars_n.u[i,j]^2 for i in 1:OCPI.meta.nu, j in 1:OCPI.settings.N-1)
reg = 0.0
obj += 0.01 * reg

# populate objective
JuMP.@objective(OCPI.model, Min, obj)

# solve NLP
JuMP.optimize!(OCPI.model)

# collect results
MTk2JuMP.IF.build.collect_all_results!(OCPI)

# collect solver info
MTk2JuMP.IF.build.collect_solver_info!(OCPI, reg=reg)

# rebuild time vector if tf is optimized
t = LinRange(OCPI.settings.tspan[1], OCPI.settings.tspan[2], OCPI.settings.N)
    
# Align the comparison plots on the same time grid used by the NLP outputs.
ty = t[1:end-1]

# Compare the values used inside the optimization against direct LUT evaluation.
# The left column shows the two signals overlaid, and the right column shows the
# pointwise difference so the interpolation drift is visible at a glance.
soc_samples = OCPI.res.x[:SOC][1:end-1]
omega_samples = OCPI.res.y[:omega_shaft]
te_samples = OCPI.res.u[:T_e]

voc_opt = OCPI.res.y[:V_oc]
r0_opt = OCPI.res.y[:R_0]
bsfc_opt = OCPI.res.y[:bsfc]

voc_ref = [OCPI.LUTs1D[:voc_lut].itp(s) for s in soc_samples]
r0_ref = [OCPI.LUTs1D[:r0_lut].itp(s) for s in soc_samples]
bsfc_ref = [OCPI.LUTs2D[:bsfc_lut].itp(ω, τ) for (ω, τ) in zip(omega_samples, te_samples)]

function comparison_panel(time, optimized, reference; ylabel::String, title::String)
    comparison = plot(
        time,
        optimized;
        label="inside optimization",
        linewidth=2,
        color=:darkblue,
        ylabel=ylabel,
        title=title,
        legend=:topright,
        grid=true,
    )
    plot!(comparison, time, reference; label="direct LUT interpolation", linewidth=2, linestyle=:dash, color=:darkorange)

    delta = plot(
        time,
        optimized .- reference;
        label="difference",
        linewidth=2,
        color=:darkred,
        ylabel="Δ" * ylabel,
        xlabel="Time \$t\$ [s]",
        legend=:topright,
        grid=true,
    )

    return comparison, delta
end

voc_cmp, voc_diff = comparison_panel(ty, voc_opt, voc_ref; ylabel="Open Circuit Voltage \$V_{oc}\$", title="\$V_{oc}\$: optimization vs LUT")
r0_cmp, r0_diff = comparison_panel(ty, r0_opt, r0_ref; ylabel="Resistance \$R_0\$", title="\$R_0\$: optimization vs LUT")
bsfc_cmp, bsfc_diff = comparison_panel(ty, bsfc_opt, bsfc_ref; ylabel="Brake Specific Fuel Consumption \$bsfc\$", title="BSFC LUT: optimization vs LUT")

# Plot output trajectories
fmt = (linewidth=2, grid=true, xlabel="Time \$t\$ [s]")

p1 = plot(t, OCPI.res.x[:x]; ylabel="Position \$x\$", legend=false, fmt...)
p2 = plot(t, OCPI.res.x[:SOC]; ylabel="State of Charge \$SOC\$", legend=false, fmt...)
p3 = plot(t, OCPI.res.x[:m_fuel]; ylabel="Fuel Mass \$m_{fuel}\$", legend=false, fmt...)
p4 = plot(t, OCPI.res.x[:v]; ylabel="Velocity \$v\$", legend=false, fmt...)
plot!(p4, t, vel_itp.(t); label="WLTP reference", linewidth=2, linestyle=:dash, color=:darkorange)

p5 = plot(t[1:end-1], OCPI.res.u[:T_e]; ylabel="Engine Torque \$T_e\$", color=:darkred, legend=false, fmt...)
p6 = plot(t[1:end-1], OCPI.res.u[:T_m]; ylabel="Motor Torque \$T_m\$", color=:darkblue, legend=false, fmt...)
p7 = plot(t[1:end-1], OCPI.res.u[:F_brk]; ylabel="Braking Force \$F_{brk}\$", color=:darkgreen, legend=false, fmt...)

plot(
    p1, p2,
    p3, p4,
    p5, p6, p7,
    voc_cmp, voc_diff,
    r0_cmp, r0_diff,
    bsfc_cmp, bsfc_diff;
    layout = (4, 4),
    size = (1920, 1080),
    dpi = 200,
    plot_title = "HEV Optimal Control and Smooth LUT Comparison",
    left_margin = 5Plots.mm,
    bottom_margin = 5Plots.mm,
)



