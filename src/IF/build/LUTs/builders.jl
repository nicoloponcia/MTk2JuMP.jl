struct SmoothRBF1D
    x::Vector{Float64}
    coeffs::Vector{Float64}
    epsilon::Float64
end

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

function _build_smooth_rbf1d(x::Vector{Float64}, y::Vector{Float64}; shape::Union{Nothing,Float64}=nothing, reg::Float64=1e-12)
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

function build_LUT1d(x::Vector{Float64}, y::Vector{Float64}, OCPI::Any, LUTName::Symbol=:UnknownLUT)
    itp = _build_smooth_rbf1d(x, y)
    itp_wrapper(t) = itp(t)
    first_der(t) = ForwardDiff.derivative(itp_wrapper, t)
    second_der(t) = ForwardDiff.derivative(first_der, t)

    itp_operator = add_nonlinear_operator(
        OCPI.model,
        1,
        itp_wrapper,
        first_der,
        second_der;
        name=LUTName,
    )

    return LUT1D(x, y, itp, itp_wrapper, itp_operator)
end

function build_LUT2d(x1::Vector{Float64}, x2::Vector{Float64}, y::Matrix{Float64}, OCPI::OCPInterface_, LUTName::Symbol)
    itp2d = DataInterpolations2D.build_smooth_rbf2d(x1, x2, y)
    itp2d_wrapper(t1, t2) = itp2d(t1, t2)
    f_vec(v) = itp2d_wrapper(v[1], v[2])

    function first_der(g::AbstractVector, t1, t2)
        v = SVector(t1, t2)
        g .= ForwardDiff.gradient(f_vec, v)
        return nothing
    end

    function second_der(H::AbstractMatrix, t1, t2)
        v = SVector(t1, t2)
        H .= ForwardDiff.hessian(f_vec, v)
        return nothing
    end

    itp_operator = add_nonlinear_operator(
        OCPI.model,
        2,
        itp2d_wrapper,
        first_der,
        second_der;
        name=LUTName,
    )

    return LUT2D(x1, x2, y, itp2d, itp2d_wrapper, itp_operator)
end