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


tf = @variable(OCPI.model, 0 <= tf <= OCPI.settings.tspan[2] * 2, start = OCPI.settings.tspan[2])

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
@constraint(OCPI.model, OCPI.vars.x[x_idx, 1] == 100.0)
# final target on base

y_idx = findfirst(x -> x == "y(t)", OCPI.meta.x_names)
@constraint(OCPI.model, OCPI.vars.x[y_idx, 1] == 100.0)
# final target on base

vx_idx = findfirst(x -> x == "v_x(t)", OCPI.meta.x_names)
@constraint(OCPI.model, OCPI.vars.x[vx_idx, 1] == 0.0)
@constraint(OCPI.model, OCPI.vars.x[vx_idx, end] == 0.0)

vy_idx = findfirst(x -> x == "v_y(t)", OCPI.meta.x_names)
@constraint(OCPI.model, OCPI.vars.x[vy_idx, 1] == -5.0)
@constraint(OCPI.model, OCPI.vars.x[vy_idx, end] == 0.0)

θ_idx = findfirst(x -> x == "θ(t)", OCPI.meta.x_names)
@constraint(OCPI.model, OCPI.vars.x[θ_idx, 1] == pi/6)
@constraint(OCPI.model, OCPI.vars.x[θ_idx, end] == 0.0)

ω_idx = findfirst(x -> x == "ω(t)", OCPI.meta.x_names)
@constraint(OCPI.model, OCPI.vars.x[ω_idx, 1] == 0.2)
@constraint(OCPI.model, OCPI.vars.x[ω_idx, end] == 0.0)



# set auxiliary constraints
# base always above ground
y_base_expr = MTk2JuMP.IF.build.get_yi_by_name(OCPI, :y_base)
@constraint(OCPI.model, [i in 1:OCPI.settings.N-1], y_base_expr[i] >= 0.0)

# end condition on base
@constraint(OCPI.model, y_base_expr[end] == 0.0)

x_base_expr = MTk2JuMP.IF.build.get_yi_by_name(OCPI, :x_base)
@constraint(OCPI.model, x_base_expr[end] == 0.0)


# define objective to minimize energy consumption
# P_expr = MTk2JuMP.IF.build.get_yi_by_name(OCPI, :P)
# obj = sum(OCPI.vars_n.u[i,j]^2 for i in 1:OCPI.meta.nu, j in 1:OCPI.settings.N-1)
obj = tf

# add regularization on control imputs
reg = sum(OCPI.vars_n.u[i,j]^2 for i in 1:OCPI.meta.nu, j in 1:OCPI.settings.N-1)
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
t = LinRange(OCPI.settings.tspan[1], value(tf), OCPI.settings.N)

# Plot output trajectories
fmt = (linewidth=2, legend=false, grid=true, xlabel="Time \$t\$ [s]")
# Construct individual subplots with mathematical LaTeX labels
p1 = plot(t, OCPI.res.x[:x]; ylabel="Position \$x\$", fmt...)
p2 = plot(t, OCPI.res.x[:y]; ylabel="Position \$y\$", fmt...)
p3 = plot(t, OCPI.res.x[:v_x]; ylabel="Velocity \$v_x\$", fmt...)
p4 = plot(t, OCPI.res.x[:v_y]; ylabel="Velocity \$v_y\$", fmt...)
p5 = plot(t, OCPI.res.x[:θ]; ylabel="Angle \$\\theta\$", fmt...)
p6 = plot(t, OCPI.res.x[:ω]; ylabel="Ang. Vel. \$\\omega\$", fmt...)

# Control and output variables (note the distinct time vector indices)
p7 = plot(t[1:end-1], OCPI.res.u[:F_m]; ylabel="Main Thrust \$F_m\$", color=:darkred, fmt...)
p8 = plot(t[1:end-1], OCPI.res.u[:F_l]; ylabel="Lateral Thrust \$F_l\$", color=:darkgreen, fmt...)
p9 = plot(t[1:end-1], OCPI.res.u[:F_r]; ylabel="Lateral Thrust \$F_r\$", color=:darkblue, fmt...)

# Assemble all subplots into a master figure
plot(p1, p2, p3, p4, p5, p6, p7, p8, p9,
     layout = (3, 3), size = (900, 450), dpi = 200,
     plot_title = "Optimal Control Trajectories",
     left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)





# Create animation

function generate_rocket_animation(t_vec, x_vec, y_vec, θ_vec, Fm_vec, Fl_vec, Fr_vec; L=15.0, k_fm=0.0005, k_fl=0.005, skip_steps=1)
    # Determine plot limits based on data range
    x_min = minimum(x_vec) - L * 2.5
    x_max = maximum(x_vec) + L * 2.5
    y_min = -L
    y_max = maximum(y_vec) + L * 2.0

    # Enforce a minimum width to maintain a reasonable aspect ratio
    if x_max - x_min < 50.0
        center_x = (x_max + x_min) / 2.0
        x_min = center_x - 25.0
        x_max = center_x + 25.0
    end

    idxs = 1:skip_steps:length(t_vec)
    xx = x_vec[idxs]
    yy = y_vec[idxs]
    θθ = θ_vec[idxs]
    FFm = Fm_vec[idxs]
    FFl = Fl_vec[idxs]
    FFr = Fr_vec[idxs]
    tt = t_vec[idxs]

    anim = @animate for i in 1:length(idxs)
        x_c = xx[i]
        y_c = yy[i]
        θ   = θθ[i]
        F_m = FFm[i]
        F_l = FFl[i]
        F_r = FFr[i]

        # Calculate rocket geometry parameters
        sinθ = sin(θ)
        cosθ = cos(θ)
        
        x_base = x_c - L * sinθ
        y_base = y_c - L * cosθ
        
        x_nose = x_c + L * sinθ
        y_nose = y_c + L * cosθ

        # Initialize the frame
        p = plot(xlims=(x_min, x_max), ylims=(y_min, y_max), 
                 aspect_ratio=:equal, legend=false, grid=true,
                 title="Rocket Landing Trajectory (t = $(round(tt[i], digits=2))s)")
        
        # Draw the ground
        hline!(p, [0.0], color=:black, linewidth=2.0)

        # Draw the target landing pad
        plot!(p, [-10.0, 10.0], [0.0, 0.0], linewidth=6, color=:green)

        # Draw the rocket body
        plot!(p, [x_base, x_nose], [y_base, y_nose], linewidth=6, color=:gray)

        # Draw the Center of Mass
        scatter!(p, [x_c], [y_c], marker=:circle, markersize=6, color=:blue)
        
        # Draw the actuation force vectors (Exhaust plumes)
        # Main Engine Plume (Exhaust points opposite to force direction)
        if F_m > 1e-3
            quiver!(p, [x_base], [y_base], 
                    quiver=([-F_m * k_fm * sinθ], [-F_m * k_fm * cosθ]), 
                    color=:orange, linewidth=4)
        end
        
        # Left Lateral Plume (Force pushes right, exhaust points left)
        if F_l > 1e-3
            quiver!(p, [x_base], [y_base], 
                    quiver=([-F_l * k_fl * cosθ], [F_l * k_fl * sinθ]), 
                    color=:red, linewidth=2)
        end

        # Right Lateral Plume (Force pushes left, exhaust points right)
        if F_r > 1e-3
            quiver!(p, [x_base], [y_base], 
                    quiver=([F_r * k_fl * cosθ], [-F_r * k_fl * sinθ]), 
                    color=:red, linewidth=2)
        end
    end
    
    return anim
end

target_fps = 30
dt = OCPI.settings.di 
skip = max(1, round(Int, 1.0 / (target_fps * dt)))

# Generate and export the animation
# A leading zero is concatenated if the control vectors (u) are one element shorter than state vectors (x)
animation_obj = generate_rocket_animation(
    t, 
    OCPI.res.x[:x], 
    OCPI.res.x[:y], 
    OCPI.res.x[:θ], 
    vcat(0.0, OCPI.res.u[:F_m]), 
    vcat(0.0, OCPI.res.u[:F_l]), 
    vcat(0.0, OCPI.res.u[:F_r]); 
    L=15.0, 
    skip_steps=skip
)

gif(animation_obj, "examples/rocket_landing/rocket_optimal_control.gif", fps=target_fps)