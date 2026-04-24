using MTk2JuMP
using Plots

# Short alias for collocation and Lagrange polynomial utilities.
const LG = MTk2JuMP.IF.LGPoly

# Collocation order used for each method.
# All methods build a local polynomial with N + 1 support points after
# prepending the interval-start node tau = 0.
n_order = 4

# All supported Legendre-family collocation variants.
methods = [:LG, :LGR, :LGL]
method_label = Dict(
    :LG => "Legendre-Gauss (LG)",
    :LGR => "Legendre-Gauss-Radau (LGR)",
    :LGL => "Legendre-Gauss-Lobatto (LGL)",
)

# Build full transcription artifacts for each method.
systems = Dict{Symbol, NamedTuple{(:nodes, :C, :B, :D), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}}}}()
for m in methods
    nodes, C, B, D = LG.get_collocation_system(n_order, method=m)
    systems[m] = (nodes=nodes, C=C, B=B, D=D)
end

println("--- Legendre Collocation Comparison (N = $n_order) ---")

# Educational summary of what differs in direct transcription.
println("\nTranscription differences:")
println("- nodes: collocation support points in [0, 1] (with tau=0 prepended)")
println("- C: differentiation matrix used in dynamic defect equations")
println("- B: quadrature weights used in integral approximation")
println("- D: endpoint reconstruction vector used for interval closure at tau=1")

for m in methods
    sys = systems[m]
    println("\n$(method_label[m])")
    println("  nodes       = ", round.(sys.nodes; digits=8))
    println("  sum(B)      = ", round(sum(sys.B); digits=12), " (should be approximately 1)")
    println("  D[end] info = ", round(sys.D[end]; digits=8), " (larger when endpoint behavior is emphasized)")
    println("  size(C)     = ", size(sys.C), ", size(B) = ", size(sys.B), ", size(D) = ", size(sys.D))
end

println("\nMethod-specific endpoint structure:")
println("- LG  : interior Gauss points only (no forced endpoint node at tau=1)")
println("- LGR : includes right endpoint tau=1, not left endpoint (left is prepended globally)")
println("- LGL : includes both interval endpoints in the Lobatto family before assembly")

# Dense plotting grid for smooth basis-curve visualization.
tau_vec = LinRange(0.0, 1.0, 200)

# Plot 1: Basis polynomials for each method side-by-side.
basis_plots = Any[]
for m in methods
    sys = systems[m]
    nodes = sys.nodes

    p_basis = plot(
        xlabel = "tau (-)",
        ylabel = "L_i(tau)",
        title = "$(method_label[m]) Basis (N = $n_order)",
        legend = :outertop,
        legend_columns = 3,
        dpi = 110,
        xlim = (0, 1),
        size = (780, 360),
    )

    hline!(p_basis, [0.0], linestyle=:dash, color=:black, alpha=0.5, linewidth=1.0, label="")
    hline!(p_basis, [1.0], linestyle=:dash, color=:red, alpha=0.5, linewidth=1.0, label="")

    colors = palette(:viridis, length(nodes))
    for i in eachindex(nodes)
        L_vals = [LG.evaluate_lagrange_polynomial(i, tau, nodes) for tau in tau_vec]
        plot!(p_basis, tau_vec, L_vals, color=colors[i], label="L_$(i)")
        vline!(p_basis, [nodes[i]], linestyle=:dot, color=colors[i], alpha=0.35, linewidth=1.2, label="")
        scatter!(p_basis, [nodes[i]], [1.0], color=colors[i], markerstrokecolor=:black, markersize=4, label="")
    end

    push!(basis_plots, p_basis)
end

display(plot(basis_plots..., layout=(1, 3), size=(1600, 430), dpi=120))

# Plot 2: Direct transcription comparison of nodes, B, and D.
p_nodes = plot(
    xlabel="tau (-)",
    ylabel="method index",
    title="Collocation Nodes by Method",
    yticks=(1:3, [string(m) for m in methods]),
    xlim=(0, 1),
    legend=false,
    dpi=120,
)

for (k, m) in enumerate(methods)
    x = systems[m].nodes
    y = fill(k, length(x))
    scatter!(p_nodes, x, y, markersize=7)
end

p_B = plot(
    xlabel="basis index i",
    ylabel="B_i",
    title="Quadrature Weights (B)",
    legend=:topright,
    dpi=120,
)
for m in methods
    plot!(p_B, 1:length(systems[m].B), systems[m].B, marker=:circle, linewidth=2.0, label=string(m))
end

p_D = plot(
    xlabel="basis index i",
    ylabel="D_i",
    title="Endpoint Reconstruction (D)",
    legend=:topright,
    dpi=120,
)
for m in methods
    plot!(p_D, 1:length(systems[m].D), systems[m].D, marker=:diamond, linewidth=2.0, label=string(m))
end

display(plot(p_nodes, p_B, p_D, layout=(1, 3), size=(1600, 430), dpi=120))

# Print C matrices for exact transcription inspection.
println("\nDetailed differentiation matrices C:")
for m in methods
    println("\n$(method_label[m]) -> C")
    display(round.(systems[m].C; digits=8))
end
