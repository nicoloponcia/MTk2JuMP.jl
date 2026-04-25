struct SmoothRBF1D
    x::Vector{Float64}
    coeffs::Vector{Float64}
    epsilon::Float64
end

struct SmoothRBF2D
    x1::Vector{Float64}
    x2::Vector{Float64}
    coeffs::Vector{Float64}
    px::Vector{Float64}
    py::Vector{Float64}
    epsilon::Float64
end
