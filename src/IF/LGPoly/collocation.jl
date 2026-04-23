function compute_differentiation_vector(nodes::Vector{Float64})
    N = length(nodes)
    D = zeros(N)

    tau_final = 1.0

    for j in 1:N
        lagrange_val = 1.0
        for k in 1:N
            if k != j
                lagrange_val *= (tau_final - nodes[k]) / (nodes[j] - nodes[k])
            end
        end
        D[j] = lagrange_val
    end

    return D
end

function get_collocation_system(N::Int; method::Symbol=:LGR)
    N < 1 && throw(ArgumentError("Number of collocation nodes N must be >= 1, got N=$N"))
    method ∉ (:LG, :LGR, :LGL) && throw(ArgumentError("Unknown method $method. Use :LG, :LGR, or :LGL"))

    collocation_nodes = if method == :LG
        get_lg_nodes(N)
    elseif method == :LGR
        get_lgr_nodes(N)
    else
        get_lgl_nodes(N)
    end

    nodes = [0.0; collocation_nodes]

    C = compute_differentiation_matrix(nodes)
    B = compute_integration_weights(nodes)
    D = compute_differentiation_vector(nodes)

    return nodes, C, B, D
end

function compute_differentiation_matrix(nodes::Vector{Float64})
    N = length(nodes)
    C = zeros(N, N)

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

    for i in 1:N
        for j in 1:N
            if i != j
                C[i, j] = (weights[j] / weights[i]) / (nodes[i] - nodes[j])
            end
        end
        C[i, i] = -sum(C[i, k] for k in 1:N if k != i)
    end

    return C
end

function compute_integration_weights(nodes::Vector{Float64})
    N = length(nodes)
    B = zeros(N)

    gl_order = N + 5
    gl_nodes, gl_weights = gauss_legendre_quadrature(gl_order)

    for j in 1:N
        integral = 0.0
        for k in 1:gl_order
            tau = 0.5 * (gl_nodes[k] + 1.0)
            weight = 0.5 * gl_weights[k]

            lagrange_val = evaluate_lagrange_polynomial(j, tau, nodes)
            integral += weight * lagrange_val
        end
        B[j] = integral
    end

    return B
end

function gauss_legendre_quadrature(n::Int)
    n < 1 && throw(ArgumentError("Number of quadrature points n must be >= 1, got n=$n"))

    nodes = zeros(n)
    weights = zeros(n)

    for i in 1:n
        x = cos(pi * (i - 0.25) / (n + 0.5))

        for _ in 1:10
            P, dP = legendre_poly(n, x)
            x -= P / dP
        end

        nodes[i] = x

        _, dP = legendre_poly(n, x)
        weights[i] = 2.0 / ((1 - x^2) * dP^2)
    end

    return nodes, weights
end
