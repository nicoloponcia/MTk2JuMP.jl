# module types

mutable struct CollSettings_
    poly::Symbol
    order::Int
    nodes::Vector{Float64}
    B::Vector{Float64}
    C::Matrix{Float64}
    D::Vector{Float64}
end
CollSettings(; poly::Symbol=:LGR, order::Int=3) = CollSettings_(poly, order, Float64[], Float64[], zeros(Float64,0,0), Float64[])

mutable struct OCPSettings_
    N::Int
    tspan::Tuple{Float64,Float64}
    di::Float64
    int_method::Symbol
    Coll_set::CollSettings_
    Params_mode::Symbol
end
OCPSettings() = OCPSettings_(0, (0.0,0.0), 0.0, :None, CollSettings(), :None)

mutable struct gDyn_
    g::Array{JuMP.ConstraintRef}
    lhs::Array{JuMP.NonlinearExpr}
    rhs::Array{JuMP.NonlinearExpr}
    closure::Array{JuMP.ConstraintRef}
end
gDyn() = gDyn_(JuMP.ConstraintRef[], JuMP.NonlinearExpr[], JuMP.NonlinearExpr[], JuMP.ConstraintRef[])


mutable struct metadata_
    nx::Int
    nu::Int
    np::Int
    ny::Int

    x_names ::Vector{String}
    u_names ::Vector{String}
    p_names ::Vector{String} 
    y_names ::Vector{String}
end
metadata() = metadata_(0, 0, 0, 0, String[], String[], String[], String[])


mutable struct vars_n_
    x::Array{JuMP.VariableRef}
    u::Array{JuMP.VariableRef}
    p::Array{JuMP.VariableRef}
    
    x_col::Array{JuMP.VariableRef}
    u_col::Array{JuMP.VariableRef}
end
vars_n() = vars_n_(JuMP.VariableRef[], JuMP.VariableRef[], JuMP.VariableRef[], JuMP.VariableRef[], JuMP.VariableRef[])

mutable struct vars_
    x::Array{JuMP.AffExpr}
    u::Array{JuMP.AffExpr}
    p::Array{JuMP.AffExpr}

    x_col::Array{JuMP.AffExpr}
    u_col::Array{JuMP.AffExpr}
end
vars() = vars_(JuMP.AffExpr[], JuMP.AffExpr[], JuMP.AffExpr[], JuMP.AffExpr[], JuMP.AffExpr[])

# mutable struct vars_s_
#     x::Array{Symbolics.Num}
#     u::Array{Symbolics.Num}
#     p::Array{Symbolics.Num}
# end
# vars_s() = vars_s_(Symbolics.Num[], Symbolics.Num[], Symbolics.Num[])


mutable struct bounds_
    xL::Array{Float64}
    xU::Array{Float64}

    uL::Array{Float64}
    uU::Array{Float64}

    pL::Array{Float64}
    pU::Array{Float64}
end
bounds() = bounds_(Float64[], Float64[], Float64[], Float64[], Float64[], Float64[])

mutable struct scales_
    x::Array{Float64}
    u::Array{Float64}
    p::Array{Float64}
    dx::Array{Float64}
end
scales() = scales_(Float64[], Float64[], Float64[], Float64[])

mutable struct res_
    x::Dict{Symbol, Vector{Float64}}
    x_col::Dict{Symbol, Matrix{Float64}}
    u::Dict{Symbol, Vector{Float64}}
    y::Dict{Symbol, Vector{Float64}}
    p::Dict{Symbol, Vector{Float64}}
end
res() = res_(Dict{Symbol, Vector{Float64}}(), Dict{Symbol, Matrix{Float64}}(), Dict{Symbol, Vector{Float64}}(),
             Dict{Symbol, Vector{Float64}}(), Dict{Symbol, Vector{Float64}}())

mutable struct solver_info_
    objective_value::Float64
    regularization_value::Float64
    solver_time::Float64
    termination_status::Any
    iterations::Int
    primals::Dict{Symbol, Vector{Float64}}
end
solver_info() = solver_info_(0.0, 0.0, 0.0, nothing, 0, Dict{Symbol, Vector{Float64}}())

mutable struct duals_
    gDyn::Dict{Symbol, Matrix{Float64}}
    gAux::Dict{Symbol, Array{Float64}}
end
duals() = duals_(Dict{Symbol, Matrix{Float64}}(), Dict{Symbol, Array{Float64}}())



mutable struct OCPInterface_
    model::JuMP.Model
    sys::Any

    settings::OCPSettings_

    meta::metadata_

    vars_n::vars_n_
    vars::vars_
    # vars_s::vars_s_
    decision_vars::Array{Symbolics.Num}

    bounds::bounds_
    scales::scales_

    # TODO check type
    y_exprs::Matrix{Any}

    gDyn :: gDyn_

    gAux::Dict{Symbol, Array{JuMP.ConstraintRef}}
    gX0::Array{JuMP.ConstraintRef}

    # TODO standardize
    uDer::Array{JuMP.NonlinearExpr}
    uAcc::Array{JuMP.NonlinearExpr}

    res::res_

    LUTs1D::Dict{Symbol, LUTs.LUT1D}
    LUTs2D::Dict{Symbol, LUTs.LUT2D}

    solver_info::solver_info_

    duals::duals_

    SLS::Dict{Symbol, Vector{Float64}}
end

OCPInterface() = OCPInterface_(JuMP.Model(), nothing,
                              OCPSettings(),
                              metadata(),
                              vars_n(), vars(),
                            #   vars_s(),
                              Symbolics.Num[],
                              bounds(),
                              scales(),
                              Matrix{Any}(undef, 0, 0),
                              gDyn(), Dict{Symbol, Array{JuMP.ConstraintRef}}(), JuMP.ConstraintRef[],
                              JuMP.NonlinearExpr[], JuMP.NonlinearExpr[],
                              res(),
                              Dict{Symbol, LUTs.LUT1D}(), Dict{Symbol, LUTs.LUT2D}(),
                              solver_info(),
                              duals(),
                              Dict{Symbol, Vector{Float64}}())

# end