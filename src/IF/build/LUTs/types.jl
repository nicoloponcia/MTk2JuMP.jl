mutable struct LUT1D
    x::Vector{Float64}
    y::Vector{Float64}
    itp::Any
    itp_wrapper::Any
    itp_operator::Any
end

mutable struct LUT2D
    x1::Vector{Float64}
    x2::Vector{Float64}
    y::Matrix{Float64}
    itp::Any
    itp_wrapper::Any
    itp_operator::Any
end
