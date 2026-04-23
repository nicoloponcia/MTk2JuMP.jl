function build1D(x::Vector{Float64}, y::Vector{Float64}, OCPI::Any, LUTName::Symbol=:UnknownLUT)
    ext_type = ExtrapolationType.Constant
    itp = LinearInterpolation(x, y, extrapolation_right=ext_type, extrapolation_left=ext_type)
    itp_wrapper(t) = itp(t)
    itp_operator = add_nonlinear_operator(OCPI.model, 1, itp_wrapper; name = LUTName)
    return LUT1D(x, y, itp, itp_wrapper, itp_operator)
end

function build2D(x1::Vector{Float64}, x2::Vector{Float64}, y::Matrix{Float64}, OCPI::Any, LUTName::Symbol)
    itp_dim = (
        LinearInterpolationDimension(x1),
        LinearInterpolationDimension(x2)
    )
    itp2d = NDInterpolation(y, itp_dim)
    itp2d_wrapper(t1, t2) = itp2d(t1, t2)
    itp_operator = add_nonlinear_operator(OCPI.model, 2, itp2d_wrapper; name = LUTName)
    return LUT2D(x1, x2, y, itp2d, itp2d_wrapper, itp_operator)
end
