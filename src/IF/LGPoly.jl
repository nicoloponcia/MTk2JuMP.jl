"""
    LGPoly

A Julia module for computing Legendre polynomial nodes and collocation systems
used in optimal control and numerical integration.

This module provides:
- Legendre polynomials and their derivatives
- Legendre-Gauss (LG), Legendre-Gauss-Radau (LGR), and Legendre-Gauss-Lobatto (LGL) nodes
- Collocation matrices for pseudospectral methods
- Lagrange polynomial evaluation using barycentric formulation

# References
- Trefethen, L. N. (2000). Spectral methods in MATLAB. SIAM.
- Boyd, J. P. (2001). Chebyshev and Fourier spectral methods. Dover Publications.
"""
module LGPoly

# Constants for numerical computations
const NEWTON_TOL = 1e-15        # Convergence tolerance for Newton-Raphson
const MAX_NEWTON_ITER = 100     # Maximum Newton-Raphson iterations
const NODE_TOL = 1e-15          # Tolerance for detecting exact node matches

"""
    legendre_poly(n::Int, x::Real) -> (P_n, dP_n)

Compute the n-th Legendre polynomial P_n(x) and its derivative dP_n(x) using
the three-term recurrence relation.

# Arguments
- `n::Int`: Degree of the Legendre polynomial (n ≥ 0)
- `x::Real`: Evaluation point

# Returns
- `P_n`: Value of the n-th Legendre polynomial at x
- `dP_n`: Value of the derivative of the n-th Legendre polynomial at x

# Mathematical Formula
The Legendre polynomials satisfy:
- P_0(x) = 1, P_1(x) = x
- (k+1)P_{k+1}(x) = (2k+1)xP_k(x) - kP_{k-1}(x)

# Example
```julia
P, dP = legendre_poly(3, 0.5)  # P_3(0.5) and P'_3(0.5)
```
"""
function legendre_poly(n::Int, x::Real)
    n < 0 && throw(ArgumentError("Polynomial degree n must be non-negative, got n=$n"))
    
    # Base cases
    if n == 0
        return 1.0, 0.0
    elseif n == 1
        return float(x), 1.0
    end
    
    # Initialize recurrence
    P_prev, P_curr = 1.0, float(x)
    dP_prev, dP_curr = 0.0, 1.0
    
    # Three-term recurrence for P_k and dP_k
    for k in 1:(n-1)
        # (k+1)P_{k+1} = (2k+1)xP_k - kP_{k-1}
        P_next = ((2k + 1) * x * P_curr - k * P_prev) / (k + 1)
        
        # Derivative recurrence: dP_{k+1} = ((2k+1)(P_k + x*dP_k) - k*dP_{k-1}) / (k+1)
        dP_next = ((2k + 1) * (P_curr + x * dP_curr) - k * dP_prev) / (k + 1)
        
        # Update for next iteration
        P_prev, P_curr = P_curr, P_next
        dP_prev, dP_curr = dP_curr, dP_next
    end
    
    return P_curr, dP_curr
end

"""
    newton_raphson_root(f, df, x0::Real; tol=NEWTON_TOL, max_iter=MAX_NEWTON_ITER) -> Real

Find a root of function f using Newton-Raphson method with derivative df.

# Arguments
- `f`: Function to find root of
- `df`: Derivative of function f  
- `x0::Real`: Initial guess
- `tol`: Convergence tolerance (default: NEWTON_TOL)
- `max_iter`: Maximum iterations (default: MAX_NEWTON_ITER)

# Returns
- Root of function f near x0

# Note
This is a helper function used internally for finding polynomial roots.
"""
function newton_raphson_root(f, df, x0::Real; tol=NEWTON_TOL, max_iter=MAX_NEWTON_ITER)
    x = float(x0)
    for iter in 1:max_iter
        fx = f(x)
        dfx = df(x)
        abs(dfx) < eps() && break  # Avoid division by zero
        
        dx = fx / dfx
        x -= dx
        
        if abs(dx) < tol
            break
        end
    end
    return x
end

"""
    get_lg_nodes(N::Int) -> Vector{Float64}

Compute the N Legendre-Gauss quadrature nodes in the interval [0,1].

These are the roots of the N-th Legendre polynomial P_N(x) transformed from [-1,1] to [0,1].

# Arguments
- `N::Int`: Number of nodes (N ≥ 1)

# Returns
- `Vector{Float64}`: N quadrature nodes in [0,1], sorted in ascending order

# Mathematical Background
Legendre-Gauss nodes are optimal for polynomial interpolation and quadrature.
They are the roots of P_N(x) and provide exact integration for polynomials
of degree up to 2N-1.

# Example
```julia
nodes = get_lg_nodes(5)  # 5 LG nodes in [0,1]
```
"""
function get_lg_nodes(N::Int)
    N < 1 && throw(ArgumentError("Number of nodes N must be ≥ 1, got N=$N"))
    
    nodes_std = zeros(N)  # Nodes on [-1,1]
    
    # Find roots of P_N(x) using Newton-Raphson with Chebyshev initial guesses
    for i in 1:N
        # Chebyshev-Gauss approximation for initial guess
        x_init = -cos(π * (i - 0.25) / (N + 0.5))
        
        # Define functions for Newton-Raphson
        f(x) = legendre_poly(N, x)[1]      # P_N(x)
        df(x) = legendre_poly(N, x)[2]     # P'_N(x)
        
        # Find root
        root = newton_raphson_root(f, df, x_init)
        nodes_std[i] = root
    end
    
    # Transform from [-1,1] to [0,1]: x_new = (x_old + 1) / 2
    return (nodes_std .+ 1.0) ./ 2.0
end

"""
    get_lgr_nodes(N::Int) -> Vector{Float64}

Compute the N Legendre-Gauss-Radau quadrature nodes in the interval [0,1].

These nodes include the right endpoint (τ=1) and are roots of P_N(x) - P_{N-1}(x).
LGR methods are commonly used in optimal control for their excellent stability properties.

# Arguments  
- `N::Int`: Number of collocation nodes (N ≥ 1)

# Returns
- `Vector{Float64}`: N quadrature nodes in [0,1], with the last node being exactly 1.0

# Mathematical Background
LGR nodes are optimal for problems where boundary conditions are specified 
at the final time. They provide excellent convergence for direct collocation methods.

# Example
```julia
nodes = get_lgr_nodes(4)  # 4 LGR nodes, last one is 1.0
```
"""
function get_lgr_nodes(N::Int)
    N < 1 && throw(ArgumentError("Number of nodes N must be ≥ 1, got N=$N"))
    
    nodes_std = zeros(N)  # Nodes on [-1,1] 
    nodes_std[end] = 1.0  # Right endpoint is included
    
    # Initial guesses using modified Chebyshev distribution
    for i in 1:(N-1)
        nodes_std[i] = -cos(2π * i / (2N + 1))
    end
    
    # Find roots of (P_N - P_{N-1})(x) using Newton-Raphson
    for i in 1:(N-1)
        # Define function and derivative for Newton-Raphson
        f(x) = begin
            pn, _ = legendre_poly(N, x)
            pnm1, _ = legendre_poly(N-1, x) 
            return pn - pnm1
        end
        
        df(x) = begin
            _, dpn = legendre_poly(N, x)
            _, dpnm1 = legendre_poly(N-1, x)
            return dpn - dpnm1
        end
        
        # Find root
        root = newton_raphson_root(f, df, nodes_std[i])
        nodes_std[i] = root
    end
    
    sort!(nodes_std)
    
    # Transform from [-1,1] to [0,1]: x_new = (x_old + 1) / 2
    return (nodes_std .+ 1.0) ./ 2.0
end

"""
    get_lgl_nodes(N::Int) -> Vector{Float64}

Compute the N Legendre-Gauss-Lobatto quadrature nodes in the interval [0,1].

These nodes include both endpoints (τ=0 and τ=1) and are roots of the derivative
of the (N-1)-th Legendre polynomial P'_{N-1}(x). LGL methods are useful when
boundary conditions are specified at both initial and final times.

# Arguments
- `N::Int`: Number of nodes (N ≥ 2, includes both endpoints)

# Returns  
- `Vector{Float64}`: N quadrature nodes in [0,1], with first node being 0.0 and last being 1.0

# Mathematical Background
LGL nodes are optimal for problems with boundary conditions at both endpoints.
The interior nodes are roots of P'_{N-1}(x), where the derivative can be computed
using the Legendre differential equation.

# Example
```julia  
nodes = get_lgl_nodes(5)  # 5 LGL nodes: 0.0, ..., ..., ..., 1.0
```
"""
function get_lgl_nodes(N::Int)
    N < 2 && throw(ArgumentError("LGL nodes require N ≥ 2 (need both endpoints), got N=$N"))
    
    nodes_std = zeros(N)  # Nodes on [-1,1]
    
    # Endpoints are fixed  
    nodes_std[1] = -1.0
    nodes_std[end] = 1.0
    
    # Base case: only two nodes
    if N == 2
        return (nodes_std .+ 1.0) ./ 2.0
    end
    
    # Find interior roots of P'_{N-1}(x) using Newton-Raphson
    poly_degree = N - 1
    
    for i in 2:(N-1)
        # Chebyshev-Gauss-Lobatto initial guess
        x_init = -cos(π * (i - 1) / (N - 1))
        
        # Define derivative function and its derivative for Newton-Raphson
        f(x) = legendre_poly(poly_degree, x)[2]  # P'_{N-1}(x)
        
        df(x) = begin
            P, dP = legendre_poly(poly_degree, x)
            # Use Legendre differential equation: (1-x²)P'' - 2xP' + n(n+1)P = 0
            # => P'' = (2xP' - n(n+1)P) / (1 - x²)
            return (2x * dP - poly_degree * (poly_degree + 1) * P) / (1 - x^2)
        end
        
        # Find root
        root = newton_raphson_root(f, df, x_init)
        nodes_std[i] = root
    end
    
    # Transform from [-1,1] to [0,1]: x_new = (x_old + 1) / 2  
    return (nodes_std .+ 1.0) ./ 2.0
end

# ============================================================================
# COLLOCATION SYSTEM FUNCTIONS
# ============================================================================

"""
    compute_differentiation_vector(nodes::Vector{Float64}) -> Vector{Float64}

Compute the differentiation vector D for evaluating Lagrange polynomials at τ = 1.

This vector is used in pseudospectral methods to enforce continuity conditions
or evaluate the state at the final time.

# Arguments
- `nodes::Vector{Float64}`: Collocation nodes in [0,1]

# Returns
- `Vector{Float64}`: Differentiation vector D where D[j] = L_j(1)

# Mathematical Background
Each element D[j] represents the value of the j-th Lagrange basis polynomial
evaluated at the final time τ = 1.

# Example
```julia
nodes = [0.0, 0.5, 1.0]
D = compute_differentiation_vector(nodes)  # [0.0, -1.0, 2.0]
```
"""
function compute_differentiation_vector(nodes::Vector{Float64})
    N = length(nodes)
    D = zeros(N)
    
    # Evaluate every Lagrange basis polynomial at final time τ = 1.0
    τ_final = 1.0
    
    for j in 1:N
        # Lagrange polynomial: L_j(τ) = ∏_{k≠j} (τ - nodes[k]) / (nodes[j] - nodes[k])
        lagrange_val = 1.0
        for k in 1:N
            if k != j
                lagrange_val *= (τ_final - nodes[k]) / (nodes[j] - nodes[k])
            end
        end
        D[j] = lagrange_val
    end
    
    return D
end

"""
    get_collocation_system(N::Int; method::Symbol=:LGR) -> (nodes, C, B, D)

Compute a complete pseudospectral collocation system for direct methods in optimal control.

This function returns all matrices needed for transcribing continuous-time optimal 
control problems into discrete nonlinear programming problems.

# Arguments
- `N::Int`: Number of collocation nodes (N ≥ 1)  
- `method::Symbol`: Collocation method (:LG, :LGR, or :LGL)

# Returns
A tuple containing:
- `nodes::Vector{Float64}`: Augmented nodes [0; collocation_nodes] in [0,1]
- `C::Matrix{Float64}`: Differentiation matrix (N+1 × N+1)
- `B::Vector{Float64}`: Integration weights (N+1,)  
- `D::Vector{Float64}`: Final state evaluation vector (N+1,)

# Mathematical Background
- **C matrix**: Approximates derivatives using Lagrange polynomial differentiation
- **B vector**: Provides integration weights for defect constraints  
- **D vector**: Evaluates polynomials at final time for continuity

# Methods
- `:LG`: Legendre-Gauss (interior nodes only)
- `:LGR`: Legendre-Gauss-Radau (includes final point)
- `:LGL`: Legendre-Gauss-Lobatto (includes both endpoints)

# Example
```julia
nodes, C, B, D = get_collocation_system(4, method=:LGR)
```
"""
function get_collocation_system(N::Int; method::Symbol=:LGR)
    # Validate input
    N < 1 && throw(ArgumentError("Number of collocation nodes N must be ≥ 1, got N=$N"))
    method ∉ (:LG, :LGR, :LGL) && throw(ArgumentError("Unknown method $method. Use :LG, :LGR, or :LGL"))
    
    # 1. Generate collocation nodes
    collocation_nodes = if method == :LG
        get_lg_nodes(N)
    elseif method == :LGR  
        get_lgr_nodes(N)
    else  # :LGL
        get_lgl_nodes(N)
    end
    
    # Augment with initial time τ = 0
    nodes = [0.0; collocation_nodes] 
    n_total = length(nodes)
    
    # 2. Differentiation matrix C
    C = compute_differentiation_matrix(nodes)
    
    # 3. Integration weights vector B  
    B = compute_integration_weights(nodes)

    # 4. Final evaluation vector D
    D = compute_differentiation_vector(nodes)

    return nodes, C, B, D
end

"""
    compute_differentiation_matrix(nodes::Vector{Float64}) -> Matrix{Float64}

Compute the differentiation matrix C for Lagrange polynomial collocation.

The differentiation matrix relates function values at nodes to derivative values:
dx/dt ≈ C * x, where x contains function values at the collocation nodes.

# Arguments
- `nodes::Vector{Float64}`: Collocation nodes in [0,1]

# Returns  
- `Matrix{Float64}`: Differentiation matrix C (N×N)

# Mathematical Background
Uses the barycentric formulation for numerical stability:
C[i,j] = (w_j/w_i) / (nodes[i] - nodes[j]) for i≠j
C[i,i] = -∑(C[i,k]) for k≠i
where w_i are the barycentric weights.
"""
function compute_differentiation_matrix(nodes::Vector{Float64})
    N = length(nodes)
    C = zeros(N, N)
    
    # Compute barycentric weights
    weights = zeros(N)
    for i in 1:N
        weight_val = 1.0
        for k in 1:N
            if k != i
                weight_val *= (nodes[i] - nodes[k])
            end
        end
        weights[i] = 1.0 / weight_val
    end
    
    # Fill off-diagonal elements using barycentric formula
    for i in 1:N
        for j in 1:N
            if i != j
                C[i, j] = (weights[j] / weights[i]) / (nodes[i] - nodes[j])
            end
        end
        # Diagonal element: negative sum of off-diagonal elements in row
        C[i, i] = -sum(C[i, k] for k in 1:N if k != i)
    end
    
    return C
end

"""
    compute_integration_weights(nodes::Vector{Float64}) -> Vector{Float64}

Compute integration weights for Lagrange polynomial quadrature.

The integration weights allow approximation of definite integrals:
∫₀¹ f(τ)dτ ≈ ∑ᵢ B[i] * f(nodes[i])

# Arguments  
- `nodes::Vector{Float64}`: Collocation nodes in [0,1]

# Returns
- `Vector{Float64}`: Integration weights B

# Mathematical Background
Each weight B[j] equals the integral of the j-th Lagrange basis polynomial:
B[j] = ∫₀¹ L_j(τ)dτ
computed using high-order Gauss-Legendre quadrature for accuracy.
"""
function compute_integration_weights(nodes::Vector{Float64})
    N = length(nodes)
    B = zeros(N)
    
    # Use high-order Gauss-Legendre quadrature for accurate integration
    # Each weight B[j] = ∫₀¹ L_j(τ)dτ where L_j is the j-th Lagrange polynomial
    gl_order = N + 5  # Extra precision
    gl_nodes, gl_weights = gauss_legendre_quadrature(gl_order) 
    
    for j in 1:N
        integral = 0.0
        for k in 1:gl_order
            # Map Gauss-Legendre node from [-1,1] to [0,1]
            τ = 0.5 * (gl_nodes[k] + 1.0)
            weight = 0.5 * gl_weights[k]
            
            # Evaluate j-th Lagrange polynomial at τ
            lagrange_val = evaluate_lagrange_polynomial(j, τ, nodes)
            integral += weight * lagrange_val
        end
        B[j] = integral
    end

    return B
end

# ============================================================================
# GAUSS-LEGENDRE QUADRATURE
# ============================================================================

"""
    gauss_legendre_quadrature(n::Int) -> (nodes, weights)

Compute nodes and weights for n-point Gauss-Legendre quadrature on [-1,1].

# Arguments
- `n::Int`: Number of quadrature points

# Returns
- `nodes::Vector{Float64}`: Quadrature nodes
- `weights::Vector{Float64}`: Quadrature weights

# Mathematical Background
Gauss-Legendre quadrature provides exact integration for polynomials of degree ≤ 2n-1.
"""
function gauss_legendre_quadrature(n::Int)
    n < 1 && throw(ArgumentError("Number of quadrature points n must be ≥ 1, got n=$n"))
    
    nodes = zeros(n)
    weights = zeros(n)
    
    # Find roots of P_n(x) and compute weights
    for i in 1:n
        # Chebyshev initial guess
        x = cos(π * (i - 0.25) / (n + 0.5))
        
        # Newton-Raphson to find root
        for iter in 1:10
            P, dP = legendre_poly(n, x)
            x -= P / dP
        end
        
        nodes[i] = x
        
        # Compute weight using formula: w_i = 2 / ((1-x_i²)[P'_n(x_i)]²)
        _, dP = legendre_poly(n, x)
        weights[i] = 2.0 / ((1 - x^2) * dP^2)
    end
    
    return nodes, weights
end

# ============================================================================
# LAGRANGE POLYNOMIAL EVALUATION
# ============================================================================

"""
    evaluate_lagrange_polynomial(i::Int, τ::Float64, nodes::Vector{Float64}) -> Float64

Evaluate the i-th Lagrange basis polynomial at point τ using barycentric formulation.

# Arguments
- `i::Int`: Index of the Lagrange polynomial (1 ≤ i ≤ length(nodes))
- `τ::Float64`: Evaluation point  
- `nodes::Vector{Float64}`: Collocation nodes

# Returns
- `Float64`: Value of L_i(τ)

# Mathematical Background  
Uses the numerically stable barycentric formulation:
L_i(τ) = [w_i/(τ-nodes[i])] / [∑_k w_k/(τ-nodes[k])]
where w_k are barycentric weights. For τ = nodes[j], returns δ_{ij}.

# Example
```julia
nodes = [0.0, 0.5, 1.0]  
val = evaluate_lagrange_polynomial(2, 0.25, nodes)  # L_2(0.25)
```
"""
function evaluate_lagrange_polynomial(i::Int, τ::Float64, nodes::Vector{Float64})
    N = length(nodes)
    
    # Validate input
    (i < 1 || i > N) && throw(BoundsError(nodes, i))
    
    # Handle exact node matches (removable singularity)
    for k in 1:N
        if abs(τ - nodes[k]) < NODE_TOL
            return k == i ? 1.0 : 0.0
        end
    end

    # Compute barycentric weights for all nodes
    weights = zeros(N)
    for j in 1:N
        weight_val = 1.0
        for k in 1:N
            if k != j
                weight_val *= (nodes[j] - nodes[k])
            end
        end
        weights[j] = 1.0 / weight_val
    end

    # Apply barycentric interpolation formula
    # L_i(τ) = [w_i/(τ-nodes[i])] / [∑_k w_k/(τ-nodes[k])]
    numerator = weights[i] / (τ - nodes[i])
    
    denominator = 0.0
    for k in 1:N
        denominator += weights[k] / (τ - nodes[k])
    end

    return numerator / denominator
end


function set_coll_method!(Coll_set)
    if Coll_set.poly ∉ (:LG, :LGR, :LGL)
        throw(ArgumentError("Unknown method $(Coll_set.poly). Use :LG, :LGR, or :LGL"))
    end

    nodes, C, B, D = get_collocation_system(Coll_set.order, method=Coll_set.poly)
    Coll_set.nodes = nodes
    Coll_set.C = C
    Coll_set.B = B
    Coll_set.D = D
end







#### TST
# using LinearAlgebra

# # include("src/MTk2JuMP.jl")
# # --- Example Usage ---

# # Define Order
# n_order = 3

# println("--- Computing LGR for N = $n_order ---")

# # 1. Get Matrices
# nodes, C, B, D = get_collocation_system(n_order, method=:LGR)

# println("\n1. Nodes (taus):")
# display(nodes)

# println("\n2. Quadrature Weights (B Vector):")
# display(B)
# println("   Sum of weights (should be 1): ", sum(B))
# println("\n3. Derivative Matrix (C Matrix):")
# display(C)
# println("\n4. Continuity Coefficients (D Vector):")
# display(D)

# plot poly
# using Plots
# using LaTeXStrings
# tau_vec = LinRange(0.0, 1.0, 200)

# p = plot(
#     xlabel = L"\tau" * " (-)", 
#     ylabel = L"L_i(\tau)", 
#     title = "Lagrange Gauss Radau Basis Polynomials (N = $n_order)",
#     legend = :outertop,
#     legend_columns = 4,
#     dpi = 100,
#     xlim = (0, 1),
#     size = (750, 600),
# )

# hline!(p, [0.0], linestyle=:dash, color=:black, alpha=0.5, linewidth=1.2, label="")
# hline!(p, [1.0], linestyle=:dash, color=:red, alpha=0.5, linewidth=1.2, label="")

# colors = palette(:viridis, n_order + 1)

# # --- 2. Tangent Line Setup ---
# target_idx = 3           # We want the neighborhood around the L3 node
# x0 = nodes[target_idx]   # The exact x-coordinate
# delta_x = 0.06           # How wide the tangent lines should be
# x_tangent = [x0 - delta_x, x0 + delta_x] # X-coordinates for the short line segments

# for i in 1:n_order+1
#     # Evaluate polynomial
#     L_vals = [evaluate_lagrange_polynomial(i, tau, nodes) for tau in tau_vec]
    
#     # Plot the polynomial curve
#     plot!(p, tau_vec, L_vals, color = colors[i], label = L"L_{%$i}(\tau)")
          
#     # Add vertical line at the specific node
#     vline!(p, [nodes[i]], linestyle = :dash, color = colors[i], alpha = 0.4, linewidth = 1.5, label = "")
           
#     # Add a marker at the peak
#     scatter!(p, [nodes[i]], [1.0], color = colors[i], markerstrokecolor = :black, markersize = 5, label = "")
    
#     # --- Calculate and Plot the Tangent Line at x0 ---
#     # 1. Find the slope using a central finite difference
#     eps_fd = 1e-5
#     y_plus  = evaluate_lagrange_polynomial(i, x0 + eps_fd, nodes)
#     y_minus = evaluate_lagrange_polynomial(i, x0 - eps_fd, nodes)
#     slope = (y_plus - y_minus) / (2 * eps_fd)
    
#     # 2. Find the exact y-value at x0
#     y0 = evaluate_lagrange_polynomial(i, x0, nodes)
    
#     # 3. Create the Y-coordinates for the tangent line using point-slope form: y = y0 + m*(x - x0)
#     y_tangent = y0 .+ slope .* (x_tangent .- x0)
    
#     # 4. Plot the short tangent segment (made slightly thicker and solid to stand out)
#     plot!(p, x_tangent, y_tangent, 
#           color = colors[i], 
#           linewidth = 4.0,       # Thicker so it stands out from the curve
#           linestyle = :solid, 
#           label = "")
          
#     # Add a small dot at the intersection point to anchor the tangent visually
#     scatter!(p, [x0], [y0], color = colors[i], markersize = 3, markerstrokewidth=0, label="")
# end

# display(p)

end  # module LGPoly