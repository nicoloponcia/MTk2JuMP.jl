function get_lg_nodes(N::Int)
    N < 1 && throw(ArgumentError("Number of nodes N must be >= 1, got N=$N"))

    nodes_std = zeros(N)

    for i in 1:N
        x_init = -cos(pi * (i - 0.25) / (N + 0.5))

        f(x) = legendre_poly(N, x)[1]
        df(x) = legendre_poly(N, x)[2]

        root = newton_raphson_root(f, df, x_init)
        nodes_std[i] = root
    end

    return (nodes_std .+ 1.0) ./ 2.0
end

function get_lgr_nodes(N::Int)
    N < 1 && throw(ArgumentError("Number of nodes N must be >= 1, got N=$N"))

    nodes_std = zeros(N)
    nodes_std[end] = 1.0

    for i in 1:(N-1)
        nodes_std[i] = -cos(2pi * i / (2N + 1))
    end

    for i in 1:(N-1)
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

        root = newton_raphson_root(f, df, nodes_std[i])
        nodes_std[i] = root
    end

    sort!(nodes_std)
    return (nodes_std .+ 1.0) ./ 2.0
end

function get_lgl_nodes(N::Int)
    N < 2 && throw(ArgumentError("LGL nodes require N >= 2 (need both endpoints), got N=$N"))

    nodes_std = zeros(N)
    nodes_std[1] = -1.0
    nodes_std[end] = 1.0

    if N == 2
        return (nodes_std .+ 1.0) ./ 2.0
    end

    poly_degree = N - 1

    for i in 2:(N-1)
        x_init = -cos(pi * (i - 1) / (N - 1))

        f(x) = legendre_poly(poly_degree, x)[2]

        df(x) = begin
            P, dP = legendre_poly(poly_degree, x)
            return (2x * dP - poly_degree * (poly_degree + 1) * P) / (1 - x^2)
        end

        root = newton_raphson_root(f, df, x_init)
        nodes_std[i] = root
    end

    return (nodes_std .+ 1.0) ./ 2.0
end
