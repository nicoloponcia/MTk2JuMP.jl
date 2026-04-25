function build_LUT1d(x::Vector{Float64}, y::Vector{Float64}, OCPI::Any, LUTName::Symbol=:UnknownLUT)
    itp = SmoothInterpolations.build_smooth_rbf1d(x, y)
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
    itp2d = SmoothInterpolations.build_smooth_rbf2d(x1, x2, y)
    itp2d_wrapper(t1, t2) = itp2d(t1, t2)
    f_vec(v) = itp2d_wrapper(v[1], v[2])

    function first_der(g::AbstractVector, t1, t2)
        v = SVector(t1, t2)
        g .= ForwardDiff.gradient(f_vec, v)
        return nothing
    end

    function second_der(H::AbstractMatrix, t1, t2)
        v = SVector(t1, t2)
        Hfull = ForwardDiff.hessian(f_vec, v)
        @inbounds for i in 1:size(H, 1)
            for j in 1:i
                H[i, j] = Hfull[i, j]
            end
        end
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