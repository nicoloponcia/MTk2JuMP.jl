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

# Load MTk model
include("model.jl")
OCPI.sys = ActiveSuspensionModel.active_suspension()

# extract sys structure
MTk2JuMP.IF.build.get_ode_architecture!(OCPI; verbose=false)

# get and set bounds
SB = include("ScalesBounds.jl")
MTk2JuMP.IF.build.set_bounds!(OCPI, SB.Bounds)
MTk2JuMP.IF.build.set_scales!(OCPI, SB.Scales)

# create decision variables and set initial guess
x0 = zeros(Float64, OCPI.meta.nx, OCPI.settings.N)
u0 = Matrix{Float64}(undef, OCPI.meta.nu, OCPI.settings.N-1)
# u0[1, :] = SB.Bounds.uU.F .* sin.(2*pi*t[1:end-1] * 10 ./ OCPI.settings.tspan[2])
u0[1, :] .= 0.0

MTk2JuMP.IF.build.set_opt_vars!(OCPI; x0=x0, u0=u0)

# NOTE since in Config.jl Params.mode = :sparse, the variables are duplicated
println("Parameter variable matrix size: ", size(OCPI.vars.p))


# NOTE fix known input values
z_r = 0.07 .- 0.05 .* cos.(2 .* pi .* t) .- 0.02 .* cos.(20 .* pi .* t)
z_r_dot = 0.1 .* pi .* sin.(2 .* pi .* t) .+ 0.4 .* pi .* sin.(20 .* pi .* t)

z_r_idx = findfirst(u -> u == "z_r(t)", OCPI.meta.u_names)
z_r_dot_idx = findfirst(u -> u == "z_r_dot(t)", OCPI.meta.u_names)
for i in 1:OCPI.settings.N-1
    fix(OCPI.vars_n.u[z_r_idx, i], z_r[i]; force=true)
    fix(OCPI.vars_n.u[z_r_dot_idx, i], z_r_dot[i]; force=true)
end


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

# set auxiliary constraints
# e.g. below, could develop another environment to handle more complex OCP
z_idx = findfirst(x -> x == "z(t)", OCPI.meta.x_names)
@constraint(OCPI.model, OCPI.vars.x[z_idx, 1] == 0.0)

z_dot_idx = findfirst(x -> x == "z_dot(t)", OCPI.meta.x_names)
@constraint(OCPI.model, OCPI.vars.x[z_dot_idx, 1] == 0.0)



# define objective to minimize energy consumption
a_expr = MTk2JuMP.IF.build.get_yi_by_name(OCPI, :a)
obj = sum(1 * a_expr[i]^2 for i in 1:OCPI.settings.N-1) # minimize discomfort

F_idx = findfirst(u -> u == "F(t)", OCPI.meta.u_names)
reg = sum(1e-5 * OCPI.vars.u[F_idx, i]^2 for i in 1:OCPI.settings.N-1) # regularize control effort

obj += reg
JuMP.@objective(OCPI.model, Min, obj)

# solve NLP
JuMP.optimize!(OCPI.model)

# collect results
MTk2JuMP.IF.build.collect_all_results!(OCPI)

# collect solver info
MTk2JuMP.IF.build.collect_solver_info!(OCPI)


# inspect optimal parameters
println("Optimal spring stiffness k: ", OCPI.res.p[:k])
println("Optimal damping coefficient c: ", OCPI.res.p[:c])


# Plot output trajectories
fmt = (linewidth=2, legend=false, grid=true, xlabel="Time \$t\$ [s]")
# Construct individual subplots with mathematical LaTeX labels
p1 = plot(t, OCPI.res.x[:z]; ylabel="Position \$z\$", fmt...)
p2 = plot(t, OCPI.res.x[:z_dot]; ylabel="Velocity \$\\dot{z}\$", fmt...)
p3 = plot(t[1:end-1], OCPI.res.y[:a]; ylabel="Acceleration \$a\$", fmt...)

# Control and output variables (note the distinct time vector indices)
p4 = plot(t[1:end-1], OCPI.res.u[:F]; ylabel="Control Force \$F\$", fmt...)

# Assemble all subplots into a master figure
plot(p1, p2, p3, p4, 
     layout = (2, 2), size = (900, 450), dpi = 200,
     plot_title = "Optimal Control Trajectories",
     left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)


# plot constraint for scaling check
p1 = plot(t[1:end-1], value.(OCPI.gDyn.rhs[:,:,1]'); title="Scaled Dynamics RHS", fmt...)
p2 = plot(t[1:end-1], value.(OCPI.gDyn.lhs[:,:,1]'); title="Scaled Dynamics LHS", fmt...)
plot(p1, p2, layout=(2,1), size=(800,400))



# Create animation
function generate_suspension_animation(t_vec, z_vec, z_r_vec, F_vec; k_f=0.001, skip_steps=1)
    # Determine vertical plot limits based on displacement range
    z_min = min(minimum(z_vec), minimum(z_r_vec)) - 0.5
    z_max = max(maximum(z_vec), maximum(z_r_vec)) + 1.0
    
    # Horizontal limits are fixed since motion is 1D (vertical)
    x_min = -1.0
    x_max = 1.0

    idxs = 1:skip_steps:length(t_vec)-1
    zz = z_vec[idxs]
    zz_r = z_r_vec[idxs]
    FF = F_vec[idxs]
    tt = t_vec[idxs]

    anim = @animate for i in 1:length(idxs)
        z_m = zz[i]
        z_road = zz_r[i]
        F = FF[i]

        # Initialize the frame
        p = plot(xlims=(x_min, x_max), ylims=(z_min, z_max), 
                 aspect_ratio=:equal, legend=false, grid=true,
                 title="Active Suspension Trajectory (t = $(round(tt[i], digits=2))s)")
        
        # Draw the road profile
        hline!(p, [z_road], color=:black, linewidth=2, linestyle=:dash)

        # Draw the suspension linkage (spring/damper representation)
        plot!(p, [0.0, 0.0], [z_road, z_m], linewidth=3, color=:gray)

        # Draw the quarter-car mass
        scatter!(p, [0.0], [z_m], marker=:square, markersize=25, color=:blue)
        
        # Draw the actuation force vector
        if abs(F) > 1e-3
            # Quiver points vertically up or down based on the sign of F
            quiver!(p, [0.0], [z_m], quiver=([0.0], [F * k_f]), color=:red, linewidth=2)
        end
    end
    
    return anim
end

# Configuration for export
target_fps = 30
# Note: corrected OCPI.settings.di to OCPI.settings.dt
dt = OCPI.settings.di 
skip = max(1, round(Int, 1.0 / (target_fps * dt)))

# Generate and export the animation
# Assuming z_r is stored in OCPI.res.x or a similar accessible structure
animation_obj = generate_suspension_animation(
    t, 
    OCPI.res.x[:z], 
    OCPI.res.u[:z_r], 
    vcat(0.0, OCPI.res.u[:F]); 
    skip_steps=skip
)

gif(animation_obj, "examples/active_suspension_codesign/AS_animation.gif", fps=target_fps)
