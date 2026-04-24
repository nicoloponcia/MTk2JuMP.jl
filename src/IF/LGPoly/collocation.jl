@inline function _normalize_collocation_method(method::Symbol)
    return method
end

@inline function _compute_barycentric_weights(nodes::Vector{Float64})
    N = length(nodes)
    weights = Vector{Float64}(undef, N)

    @inbounds for j in 1:N
        wj = 1.0
        xj = nodes[j]
        for k in 1:N
            if k != j
                wj *= (xj - nodes[k])
            end
        end
        weights[j] = inv(wj)
    end

    return weights
end

@inline function _evaluate_lagrange_barycentric(i::Int, tau::Float64, nodes::Vector{Float64}, weights::Vector{Float64})
    N = length(nodes)
    (i < 1 || i > N) && throw(BoundsError(nodes, i))

    @inbounds for k in 1:N
        if abs(tau - nodes[k]) < NODE_TOL
            return k == i ? 1.0 : 0.0
        end
    end

    numerator = weights[i] / (tau - nodes[i])
    denominator = 0.0
    @inbounds for k in 1:N
        denominator += weights[k] / (tau - nodes[k])
    end
    return numerator / denominator
end

function compute_differentiation_vector(nodes::Vector{Float64}, weights::Vector{Float64}=_compute_barycentric_weights(nodes))
    N = length(nodes)
    D = Vector{Float64}(undef, N)
    tau_final = 1.0

    @inbounds for k in 1:N
        if abs(tau_final - nodes[k]) < NODE_TOL
            fill!(D, 0.0)
            D[k] = 1.0
            return D
        end
    end

    denominator = 0.0
    @inbounds for k in 1:N
        denominator += weights[k] / (tau_final - nodes[k])
    end

    @inbounds for j in 1:N
        D[j] = (weights[j] / (tau_final - nodes[j])) / denominator
    end

    return D
end

function get_collocation_system(N::Int; method::Symbol=:LGR)
    N < 1 && throw(ArgumentError("Number of collocation nodes N must be >= 1, got N=$N"))

    method_norm = _normalize_collocation_method(method)
    method_norm ∉ (:LG, :LGR, :LGL) && throw(ArgumentError("Unknown method $method. Use :LG, :LGR, or :LGL"))

    collocation_nodes = if method_norm == :LG
        get_lg_nodes(N)
    elseif method_norm == :LGR
        get_lgr_nodes(N)
    else
        # Build Lobatto points with both endpoints, then drop the initial endpoint
        # because the global discretization start node is prepended below.
        get_lgl_nodes(N + 1)[2:end]
    end

    nodes = vcat(0.0, collocation_nodes)
    bary_weights = _compute_barycentric_weights(nodes)

    C = compute_differentiation_matrix(nodes, bary_weights)
    B = compute_integration_weights(nodes, bary_weights)
    D = compute_differentiation_vector(nodes, bary_weights)

    return nodes, C, B, D
end

function compute_differentiation_matrix(nodes::Vector{Float64}, weights::Vector{Float64}=_compute_barycentric_weights(nodes))
    N = length(nodes)
    C = zeros(Float64, N, N)

    @inbounds for i in 1:N
        xi = nodes[i]
        wi = weights[i]
        rowsum = 0.0
        for j in 1:N
            if i != j
                cij = (weights[j] / wi) / (xi - nodes[j])
                C[i, j] = cij
                rowsum += cij
            end
        end
        C[i, i] = -rowsum
    end

    return C
end

function compute_integration_weights(nodes::Vector{Float64}, weights::Vector{Float64}=_compute_barycentric_weights(nodes))
    N = length(nodes)
    B = zeros(Float64, N)

    gl_order = N + 5
    gl_nodes, gl_weights = gauss_legendre_quadrature(gl_order)

    @inbounds for j in 1:N
        integral = 0.0
        for k in 1:gl_order
            tau = 0.5 * (gl_nodes[k] + 1.0)
            wk = 0.5 * gl_weights[k]
            integral += wk * _evaluate_lagrange_barycentric(j, tau, nodes, weights)
        end
        B[j] = integral
    end

    return B
end

function gauss_legendre_quadrature(n::Int)
    n < 1 && throw(ArgumentError("Number of quadrature points n must be >= 1, got n=$n"))

    nodes = Vector{Float64}(undef, n)
    weights = Vector{Float64}(undef, n)

    m = (n + 1) ÷ 2
    @inbounds for i in 1:m
        x = cos(pi * (i - 0.25) / (n + 0.5))

        for _ in 1:MAX_NEWTON_ITER
            P, dP = legendre_poly(n, x)
            abs(dP) < eps(Float64) && break
            dx = P / dP
            x -= dx
            abs(dx) < NEWTON_TOL && break
        end

        _, dP = legendre_poly(n, x)
        w = 2.0 / ((1.0 - x * x) * dP * dP)

        left = i
        right = n + 1 - i
        nodes[left] = -x
        nodes[right] = x
        weights[left] = w
        weights[right] = w
    end

    return nodes, weights
end
