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
MTk2JuMP.IF.build.setup_optimizer!(OCPI.model, NLPConfig)

settings = (
        N = 500,
        tspan = (0, 5.0),
        int_method = NLPConfig.Integration.int_method,
        Coll_set = NLPConfig.Integration.Collocation,
        param_mode = NLPConfig.Params.mode,
)
MTk2JuMP.IF.build.set_settings!(OCPI, settings)

# Load MTk model
include("model.jl")
OCPI.sys = CartPoleModel.cart_pole()

# extract sys structure
MTk2JuMP.IF.build.get_ode_architecture!(OCPI; verbose=false)

# get and set bounds
SB = include("ScalesBounds.jl")
MTk2JuMP.IF.build.set_bounds!(OCPI, SB.Bounds)
MTk2JuMP.IF.build.set_scales!(OCPI, SB.Scales)

# create decision variables and set initial guess
x0 = zeros(Float64, OCPI.meta.nx, OCPI.settings.N)
t = LinRange(OCPI.settings.tspan[1], OCPI.settings.tspan[2], OCPI.settings.N)
u0 = Matrix{Float64}(undef, OCPI.meta.nu, OCPI.settings.N-1)
# u0[1, :] = SB.Bounds.uU.F .* sin.(2*pi*t[1:end-1] * 10 ./ OCPI.settings.tspan[2])
u0[1, :] .= 0.0

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

# set auxiliary constraints
# e.g. below, could develop another environment to handle more complex OCP
x_idx = findfirst(x -> x == "x(t)", OCPI.meta.x_names)
@constraint(OCPI.model, OCPI.vars.x[x_idx, 1] == 0.0)
@constraint(OCPI.model, OCPI.vars.x[x_idx, end] == 0.0)

v_idx = findfirst(x -> x == "v(t)", OCPI.meta.x_names)
@constraint(OCPI.model, OCPI.vars.x[v_idx, 1] == 0.0)
@constraint(OCPI.model, OCPI.vars.x[v_idx, end] == 0.0)

θ_idx = findfirst(x -> x == "θ(t)", OCPI.meta.x_names)
@constraint(OCPI.model, OCPI.vars.x[θ_idx, 1] == 0.0)
@constraint(OCPI.model, OCPI.vars.x[θ_idx, end] == pi)

ω_idx = findfirst(x -> x == "ω(t)", OCPI.meta.x_names)
@constraint(OCPI.model, OCPI.vars.x[ω_idx, 1] == 0.0)
@constraint(OCPI.model, OCPI.vars.x[ω_idx, end] == 0.0)

a_idx = findfirst(x -> x == "a(t)", OCPI.meta.x_names)
@constraint(OCPI.model, OCPI.vars.x[a_idx, 1] == 0.0)
@constraint(OCPI.model, OCPI.vars.x[a_idx, end] == 0.0)

# define objective to minimize energy consumption
P_expr = MTk2JuMP.IF.build.get_yi_by_name(OCPI, :P)
obj = sum(P_expr[i]^2 for i in 1:OCPI.settings.N-1)

# obj = sum(OCPI.vars_n.u[i,j]^2 for i in 1:OCPI.meta.nu, j in 1:OCPI.settings.N-1)

JuMP.@objective(OCPI.model, Min, obj)

# solve NLP
JuMP.optimize!(OCPI.model)

# collect results
MTk2JuMP.IF.build.collect_all_results!(OCPI)

# collect solver info
MTk2JuMP.IF.build.collect_solver_info!(OCPI)



plot(OCPI.res.x[:x])

plot(OCPI.res.x[:v])

plot(OCPI.res.x[:θ])

plot(OCPI.res.x[:ω])

plot(OCPI.res.x[:a])

plot(OCPI.res.u[:F])

plot(OCPI.res.y[:P])



plot(value.(OCPI.gDyn.rhs[:,:,1]'))







function generate_cartpole_animation(t_vec, x_vec, θ_vec, F_vec; l=1.0, k_f=0.5, skip_steps=1)
    # Determine plot limits based on data range
    x_min = minimum(x_vec) - l - 0.5
    x_max = maximum(x_vec) + l + 0.5
    y_min = -l - 0.2
    y_max = l + 0.2

    idxs = 1:skip_steps:length(t_vec)
    xx = x_vec[idxs]
    θθ = θ_vec[idxs]
    FF = F_vec[idxs]
    tt = t_vec[idxs]

    anim = @animate for i in 1:length(idxs)
        x_c = xx[i]
        θ   = θθ[i]
        F   = FF[i]

        # Calculate pole tip coordinates (θ = 0 points down)
        x_p = x_c + l * sin(θ)
        y_p = -l * cos(θ)

        # Initialize the frame
        p = plot(xlims=(x_min, x_max), ylims=(y_min, y_max), 
                 aspect_ratio=:equal, legend=false, grid=true,
                 title="Cart-Pole Trajectory (t = $(round(tt[i], digits=2))s)")
        
        # Draw the ground track (rail)
        hline!(p, [0], color=:black, linewidth=1.5)

        # Draw the pole (line from cart to tip)
        plot!(p, [x_c, x_p], [0.0, y_p], linewidth=4, color=:orange)

        # Draw the cart
        scatter!(p, [x_c], [0.0], marker=:square, markersize=12, color=:blue)
        
        # Draw the actuation force vector
        if abs(F) > 1e-3
            # The quiver function uses relative displacements
            quiver!(p, [x_c], [0.0], quiver=([F * k_f], [0.0]), color=:red, linewidth=2)
        end
    end
    
    return anim
end

# Extract data from the Optimal Control Problem output
# (Assuming OCPI.res contains a time vector, e.g., OCPI.res.t)
t_data = LinRange(OCPI.settings.tspan[1], OCPI.settings.tspan[2], OCPI.settings.N)
x_data = OCPI.res.x[:x]
θ_data = OCPI.res.x[:θ]
F_data = vcat(0.0, OCPI.res.u[:F])
target_fps = 30
dt = OCPI.settings.di
skip = max(1, round(Int, 1.0 / (target_fps * dt)))

# Generate and export the animation
animation_obj = generate_cartpole_animation(t_data, x_data, θ_data, F_data; skip_steps=skip)  # Adjust skip_steps for faster animation if needed
gif(animation_obj, "examples/cart_pole/cart_pole_optimal_control.gif", fps=target_fps)