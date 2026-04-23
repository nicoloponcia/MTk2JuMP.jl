function evaluate_lagrange_polynomial(i::Int, tau::Float64, nodes::Vector{Float64})
    N = length(nodes)

    (i < 1 || i > N) && throw(BoundsError(nodes, i))

    for k in 1:N
        if abs(tau - nodes[k]) < NODE_TOL
            return k == i ? 1.0 : 0.0
        end
    end

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

    numerator = weights[i] / (tau - nodes[i])

    denominator = 0.0
    for k in 1:N
        denominator += weights[k] / (tau - nodes[k])
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
