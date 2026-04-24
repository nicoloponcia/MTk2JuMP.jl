@inline function _map_to_unit_interval!(nodes::Vector{Float64})
    @inbounds for i in eachindex(nodes)
        nodes[i] = 0.5 * (nodes[i] + 1.0)
    end
    return nodes
end

@inline function _newton_legendre_root(n::Int, x0::Float64)
    x = x0
    for _ in 1:MAX_NEWTON_ITER
        P, dP = legendre_poly(n, x)
        abs(dP) < eps(Float64) && break
        dx = P / dP
        x -= dx
        abs(dx) < NEWTON_TOL && break
    end
    return x
end

@inline function _newton_lgr_root(n::Int, x0::Float64)
    x = x0
    nm1 = n - 1
    for _ in 1:MAX_NEWTON_ITER
        pn, dpn = legendre_poly(n, x)
        pnm1, dpnm1 = legendre_poly(nm1, x)
        f = pn - pnm1
        df = dpn - dpnm1
        abs(df) < eps(Float64) && break
        dx = f / df
        x -= dx
        abs(dx) < NEWTON_TOL && break
    end
    return x
end

@inline function _newton_lgl_inner_root(poly_degree::Int, x0::Float64)
    x = x0
    for _ in 1:MAX_NEWTON_ITER
        P, dP = legendre_poly(poly_degree, x)
        denom = 1.0 - x * x
        abs(denom) < eps(Float64) && break
        d2P = (2.0 * x * dP - poly_degree * (poly_degree + 1) * P) / denom
        abs(d2P) < eps(Float64) && break
        dx = dP / d2P
        x -= dx
        abs(dx) < NEWTON_TOL && break
    end
    return x
end

function get_lg_nodes(N::Int)
    N < 1 && throw(ArgumentError("Number of nodes N must be >= 1, got N=$N"))

    nodes_std = Vector{Float64}(undef, N)
    @inbounds for i in 1:N
        x0 = -cos(pi * (i - 0.25) / (N + 0.5))
        nodes_std[i] = _newton_legendre_root(N, x0)
    end

    sort!(nodes_std)
    return _map_to_unit_interval!(nodes_std)
end

function get_lgr_nodes(N::Int)
    N < 1 && throw(ArgumentError("Number of nodes N must be >= 1, got N=$N"))

    nodes_std = Vector{Float64}(undef, N)
    nodes_std[end] = 1.0

    @inbounds for i in 1:(N - 1)
        x0 = -cos(2pi * i / (2N + 1))
        nodes_std[i] = _newton_lgr_root(N, x0)
    end

    sort!(nodes_std)
    return _map_to_unit_interval!(nodes_std)
end

function get_lgl_nodes(N::Int)
    N < 2 && throw(ArgumentError("LGL nodes require N >= 2 (need both endpoints), got N=$N"))

    nodes_std = Vector{Float64}(undef, N)
    nodes_std[1] = -1.0
    nodes_std[end] = 1.0

    if N > 2
        poly_degree = N - 1
        @inbounds for i in 2:(N - 1)
            x0 = -cos(pi * (i - 1) / (N - 1))
            nodes_std[i] = _newton_lgl_inner_root(poly_degree, x0)
        end
    end

    sort!(nodes_std)
    return _map_to_unit_interval!(nodes_std)
end

