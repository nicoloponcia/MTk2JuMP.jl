function legendre_poly(n::Int, x::Real)
    n < 0 && throw(ArgumentError("Polynomial degree n must be non-negative, got n=$n"))

    if n == 0
        return 1.0, 0.0
    elseif n == 1
        return float(x), 1.0
    end

    P_prev, P_curr = 1.0, float(x)
    dP_prev, dP_curr = 0.0, 1.0

    for k in 1:(n-1)
        P_next = ((2k + 1) * x * P_curr - k * P_prev) / (k + 1)
        dP_next = ((2k + 1) * (P_curr + x * dP_curr) - k * dP_prev) / (k + 1)

        P_prev, P_curr = P_curr, P_next
        dP_prev, dP_curr = dP_curr, dP_next
    end

    return P_curr, dP_curr
end

function newton_raphson_root(f, df, x0::Real; tol=NEWTON_TOL, max_iter=MAX_NEWTON_ITER)
    x = float(x0)
    for _ in 1:max_iter
        fx = f(x)
        dfx = df(x)
        abs(dfx) < eps() && break

        dx = fx / dfx
        x -= dx

        if abs(dx) < tol
            break
        end
    end
    return x
end
