using MTk2JuMP
using Plots

# Access the Legendre-Gauss-Radau utilities through a short alias.
# This keeps the code readable while still being explicit about the source module.
const LG = MTk2JuMP.IF.LGPoly

# Number of collocation nodes used by the LGR scheme.
# Increasing this order generally improves approximation accuracy, but also
# increases computational cost and matrix sizes.
n_order = 3

println("--- Computing LGR for N = $n_order ---")

# Build the full collocation system for LGR.
# Returned objects:
# - nodes: collocation grid in [0, 1], augmented with tau=0
# - C: differentiation matrix, used to approximate derivatives on the grid
# - B: quadrature weights, used to approximate integrals over [0, 1]
# - D: endpoint evaluation vector, used to reconstruct value at tau=1
nodes, C, B, D = LG.get_collocation_system(n_order, method=:LGR)

println("\n1. Nodes (taus):")
display(nodes)

println("\n2. Quadrature Weights (B Vector):")
display(B)
println("   Sum of weights (should be 1): ", sum(B))

println("\n3. Derivative Matrix (C Matrix):")
display(C)

println("\n4. Continuity Coefficients (D Vector):")
display(D)

# Dense evaluation grid used only for smooth plotting of basis polynomials.
# This is independent from the collocation grid itself.
tau_vec = LinRange(0.0, 1.0, 200)

# Create the plot canvas with publication-friendly defaults.
p = plot(
    xlabel = "tau (-)",
    ylabel = "L_i(tau)",
    title = "Lagrange Gauss Radau Basis Polynomials (N = $n_order)",
    legend = :outertop,
    legend_columns = 4,
    dpi = 100,
    xlim = (0, 1),
    size = (750, 600),
)

# Reference horizontal lines to quickly see basis polynomial sign and unit level.
hline!(p, [0.0], linestyle=:dash, color=:black, alpha=0.5, linewidth=1.2, label="")
hline!(p, [1.0], linestyle=:dash, color=:red, alpha=0.5, linewidth=1.2, label="")

# One color per basis polynomial.
colors = palette(:viridis, n_order + 1)

# Choose a node where tangent behavior is inspected.
# We cap at available nodes for safety if n_order is changed.
target_idx = min(3, n_order + 1)
x0 = nodes[target_idx]
delta_x = 0.06
x_tangent = [x0 - delta_x, x0 + delta_x]

# For each Lagrange basis polynomial L_i(tau):
# 1) evaluate on a dense grid
# 2) draw the basis curve
# 3) mark its collocation node
# 4) approximate and draw tangent segment near x0
for i in 1:(n_order + 1)
    # Evaluate L_i(tau) over tau_vec for smooth rendering.
    # Mathematically, L_i is the i-th cardinal polynomial satisfying
    # L_i(nodes[j]) = delta_ij.
    L_vals = [LG.evaluate_lagrange_polynomial(i, tau, nodes) for tau in tau_vec]

    # Plot basis function and visual anchors at its node.
    plot!(p, tau_vec, L_vals, color = colors[i], label = "L_$(i)(tau)")
    vline!(p, [nodes[i]], linestyle = :dash, color = colors[i], alpha = 0.4, linewidth = 1.5, label = "")
    scatter!(p, [nodes[i]], [1.0], color = colors[i], markerstrokecolor = :black, markersize = 5, label = "")

    # Central finite-difference approximation of dL_i/dtau at x0:
    # dL_i/dtau(x0) ~= (L_i(x0 + eps) - L_i(x0 - eps)) / (2*eps)
    # This is O(eps^2) accurate for smooth functions.
    eps_fd = 1e-5
    y_plus = LG.evaluate_lagrange_polynomial(i, x0 + eps_fd, nodes)
    y_minus = LG.evaluate_lagrange_polynomial(i, x0 - eps_fd, nodes)
    slope = (y_plus - y_minus) / (2 * eps_fd)

    # Point-slope line around x0:
    # y(tau) = y0 + slope * (tau - x0)
    y0 = LG.evaluate_lagrange_polynomial(i, x0, nodes)
    y_tangent = y0 .+ slope .* (x_tangent .- x0)

    # Overlay local tangent segment and anchor point.
    # Tangent is intentionally thick so it remains visible over the basis curve.
    plot!(p, x_tangent, y_tangent, color = colors[i], linewidth = 4.0, linestyle = :solid, label = "")
    scatter!(p, [x0], [y0], color = colors[i], markersize = 3, markerstrokewidth = 0, label = "")
end

# Display in the current Julia session (REPL/VSCode/Jupyter depending on backend).
display(p)
