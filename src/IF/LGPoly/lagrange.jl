function evaluate_lagrange_polynomial(i::Int, tau::Float64, nodes::Vector{Float64})
    weights = _compute_barycentric_weights(nodes)
    return _evaluate_lagrange_barycentric(i, tau, nodes, weights)
end

function set_coll_method!(Coll_set)
    method_norm = _normalize_collocation_method(Coll_set.poly)
    if method_norm ∉ (:LG, :LGR, :LGL)
        throw(ArgumentError("Unknown method $(Coll_set.poly). Use :LG, :LGR, or :LGL"))
    end

    nodes, C, B, D = get_collocation_system(Coll_set.order, method=method_norm)
    Coll_set.poly = method_norm
    Coll_set.nodes = nodes
    Coll_set.C = C
    Coll_set.B = B
    Coll_set.D = D
end
