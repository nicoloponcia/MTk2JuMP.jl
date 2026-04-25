@inline function _gaussian_kernel(r2, epsilon2)
    return exp(-epsilon2 * r2)
end

function _default_rbf_shape_1d(x::Vector{Float64})
    if length(x) <= 1
        return 1.0
    end
    xs = sort(x)
    dx = diff(xs)
    spacing = sum(dx) / length(dx)
    return 1.0 / max(spacing, eps(Float64))
end

function _default_rbf_shape_2d(x1::Vector{Float64}, x2::Vector{Float64})
    hx = if length(x1) > 1
        xs = sort(x1)
        dx = diff(xs)
        sum(dx) / length(dx)
    else
        1.0
    end

    hy = if length(x2) > 1
        ys = sort(x2)
        dy = diff(ys)
        sum(dy) / length(dy)
    else
        1.0
    end

    scale = max(hx, hy, eps(Float64))
    return 1.0 / scale
end

function build_smooth_rbf1d(
    x::Vector{Float64},
    y::Vector{Float64};
    shape::Union{Nothing,Float64}=nothing,
    reg::Float64=1e-12,
)
    n = length(x)
    n == length(y) || throw(ArgumentError("x and y must have the same length"))
    n > 0 || throw(ArgumentError("x and y must be non-empty"))

    epsilon = isnothing(shape) ? _default_rbf_shape_1d(x) : shape
    eps2 = epsilon * epsilon

    K = Matrix{Float64}(undef, n, n)
    @inbounds for i in 1:n
        xi = x[i]
        for j in i:n
            dx = xi - x[j]
            kij = _gaussian_kernel(dx * dx, eps2)
            K[i, j] = kij
            K[j, i] = kij
        end
    end
    @inbounds for i in 1:n
        K[i, i] += reg
    end

    coeffs = K \ y
    return SmoothRBF1D(copy(x), coeffs, epsilon)
end

function (itp::SmoothRBF1D)(t::Real)
    eps2 = itp.epsilon * itp.epsilon
    Tout = promote_type(typeof(t), Float64)
    out = zero(Tout)

    @inbounds for i in eachindex(itp.x)
        dx = t - itp.x[i]
        out += itp.coeffs[i] * _gaussian_kernel(dx * dx, eps2)
    end

    return out
end

function build_smooth_rbf2d(
    x1::Vector{Float64},
    x2::Vector{Float64},
    z::Matrix{Float64};
    shape::Union{Nothing,Float64}=nothing,
    reg::Float64=1e-12,
)
    n1 = length(x1)
    n2 = length(x2)
    size(z) == (n1, n2) || throw(ArgumentError("z must have size (length(x1), length(x2))"))
    n1 > 0 && n2 > 0 || throw(ArgumentError("x1 and x2 must be non-empty"))

    epsilon = isnothing(shape) ? _default_rbf_shape_2d(x1, x2) : shape
    eps2 = epsilon * epsilon

    n = n1 * n2
    px = Vector{Float64}(undef, n)
    py = Vector{Float64}(undef, n)
    vals = Vector{Float64}(undef, n)

    idx = 0
    @inbounds for i in 1:n1
        for j in 1:n2
            idx += 1
            px[idx] = x1[i]
            py[idx] = x2[j]
            vals[idx] = z[i, j]
        end
    end

    K = Matrix{Float64}(undef, n, n)
    @inbounds for a in 1:n
        xa = px[a]
        ya = py[a]
        for b in a:n
            dx = xa - px[b]
            dy = ya - py[b]
            kij = _gaussian_kernel(dx * dx + dy * dy, eps2)
            K[a, b] = kij
            K[b, a] = kij
        end
    end

    @inbounds for i in 1:n
        K[i, i] += reg
    end

    coeffs = K \ vals
    return SmoothRBF2D(copy(x1), copy(x2), coeffs, px, py, epsilon)
end

function (itp::SmoothRBF2D)(x::Real, y::Real)
    eps2 = itp.epsilon * itp.epsilon
    Tout = promote_type(typeof(x), typeof(y), Float64)
    out = zero(Tout)

    @inbounds for i in eachindex(itp.coeffs)
        dx = x - itp.px[i]
        dy = y - itp.py[i]
        out += itp.coeffs[i] * _gaussian_kernel(dx * dx + dy * dy, eps2)
    end

    return out
end
