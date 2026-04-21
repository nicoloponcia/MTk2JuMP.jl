module build

using ModelingToolkit, JuMP
using LinearAlgebra

include("LGPoly.jl")
include("LUTs.jl")


mutable struct CollSettings
    poly::Symbol
    order::Int
    nodes::Vector{Float64}
    B::Vector{Float64}
    C::Matrix{Float64}
    D::Vector{Float64}
end
CollSettings(; poly::Symbol=:LGR, order::Int=3) = CollSettings(poly, order, Float64[], Float64[], zeros(Float64,0,0), Float64[])

mutable struct OCPSettings
    N::Int
    tspan::Tuple{Float64,Float64}
    di::Float64
    int_method::Symbol
    Coll_set::CollSettings
    Params_mode::Symbol
end
OCPSettings() = OCPSettings(0, (0.0,0.0), 0.0, :None, CollSettings(), :None)


mutable struct gDyn
    g::Array{JuMP.ConstraintRef}
    lhs::Array{JuMP.NonlinearExpr}
    rhs::Array{JuMP.NonlinearExpr}
    closure::Array{JuMP.ConstraintRef}
end
gDyn() = gDyn(JuMP.ConstraintRef[], JuMP.NonlinearExpr[], JuMP.NonlinearExpr[], JuMP.ConstraintRef[])


mutable struct OCPInterface
    model::JuMP.Model
    # sys::ModelingToolkit.System
    sys::Any

    nx::Int
    nu::Int
    np::Int
    x_names ::Vector{String}
    u_names ::Vector{String}
    p_names ::Vector{String}
    settings::OCPSettings
    x::Array{JuMP.AffExpr}
    u::Array{JuMP.AffExpr}
    p::Array{JuMP.AffExpr}
    x_n::Array{JuMP.VariableRef}
    u_n::Array{JuMP.VariableRef}
    p_n::Array{JuMP.VariableRef}
    x_col_n::Array{JuMP.VariableRef}
    x_col::Array{JuMP.AffExpr}
    u_col_n::Array{JuMP.VariableRef}
    u_col::Array{JuMP.AffExpr}
    x_sys::Array{Symbolics.Num}
    u_sys::Array{Symbolics.Num}
    p_sys::Array{Symbolics.Num}
    decision_vars::Array{Symbolics.Num}
    xL::Array{Float64}
    xU::Array{Float64}
    uL::Array{Float64}
    uU::Array{Float64}
    pL::Array{Float64}
    pU::Array{Float64}
    xScale::Array{Float64}
    uScale::Array{Float64}
    pScale::Array{Float64}
    dxScale::Array{Float64}
    y_exprs::Matrix{Any}

    gDyn :: gDyn

    gAux::Dict{Symbol, Array{JuMP.ConstraintRef}}
    gX0::Array{JuMP.ConstraintRef}

    uAcc::Array{JuMP.NonlinearExpr}

    x_res::Dict{Symbol, Vector{Float64}}
    x_col_res::Dict{Symbol, Matrix{Float64}}
    u_res::Dict{Symbol, Vector{Float64}}
    y_res::Dict{Symbol, Vector{Float64}}
    p_res::Dict{Symbol, Vector{Float64}}

    LUTs1D::Dict{Symbol, LUTs.LUT1D}
    LUTs2D::Dict{Symbol, LUTs.LUT2D}

    objective_value::Float64
    regularization_value::Float64
    solve_time::Float64
    termination_status::Any

    gDyn_dual::Dict{Symbol, Matrix{Float64}}
    gAux_dual::Dict{Symbol, Array{Float64}}

    SLS::Dict{Symbol, Vector{Float64}}

end
OCPInterface() = OCPInterface(JuMP.Model(), nothing,
                              0,0,0,
                              String[], String[], String[],
                              OCPSettings(),
                              JuMP.AffExpr[], JuMP.AffExpr[], JuMP.AffExpr[],
                              JuMP.VariableRef[], JuMP.VariableRef[], JuMP.VariableRef[],
                              JuMP.VariableRef[], JuMP.AffExpr[], JuMP.VariableRef[], JuMP.AffExpr[],
                              Symbolics.Num[], Symbolics.Num[], Symbolics.Num[], Symbolics.Num[],
                              Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[],
                              Matrix{Any}(undef, 0, 0),
                              gDyn(), Dict{Symbol, JuMP.ConstraintRef}(), JuMP.ConstraintRef[], JuMP.NonlinearExpr[],
                              Dict{Symbol, Vector{Float64}}(), Dict{Symbol, Matrix{Float64}}(), Dict{Symbol, Vector{Float64}}(), Dict{Symbol, Vector{Float64}}(), Dict{Symbol, Vector{Float64}}(),
                              Dict{Symbol, LUTs.LUT1D}(), Dict{Symbol, LUTs.LUT2D}(),
                              0.0, 0.0, 0.0, nothing, Dict{Symbol, Matrix{Float64}}(), Dict{Symbol, Array{Float64}}(),
                              Dict{Symbol, Vector{Float64}}())

function get_ode_architecture!(OCPI::OCPInterface; verbose::Bool=false)
    t = ModelingToolkit.get_iv(OCPI.sys)

    vars   = ModelingToolkit.unknowns(OCPI.sys)
    eqns   = ModelingToolkit.equations(OCPI.sys)
    params = ModelingToolkit.parameters(OCPI.sys)
    defs   = ModelingToolkit.initial_conditions(OCPI.sys)


    inputs = filter(params) do p
        istree(p) && any(isequal(t), arguments(p))
    end

    opt_params = filter(params) do p
        is_static = !(istree(p) && any(isequal(t), arguments(p)))
        has_no_default = !haskey(defs, p)
        return is_static && has_no_default
    end

    OCPI.nx = length(vars)
    OCPI.nu = length(inputs)
    OCPI.np = length(opt_params)

    OCPI.x_names = [string(v) for v in vars]
    OCPI.u_names = [string(u) for u in inputs]
    OCPI.p_names = [string(p) for p in opt_params]

    OCPI.x_sys = vars
    OCPI.u_sys = inputs
    OCPI.p_sys = opt_params
    OCPI.decision_vars = vcat(OCPI.x_sys, OCPI.u_sys, OCPI.p_sys)

    if verbose
        # Print summary
        println("OCP Architecture Summary:")
        println("  Number of state variables (nx): ", OCPI.nx)
        println("  Number of control inputs (nu): ", OCPI.nu)
        println("  Number of optimization parameters (np): ", OCPI.np)
        # Print variable names
        println("  State variable names: ", OCPI.x_names)
        println("  Control input names: ", OCPI.u_names)
        println("  Optimization parameter names: ", OCPI.p_names)
    end

end




function set_settings!(OCPI::OCPInterface, settings::NamedTuple)
    OCPI.settings.N = settings.N
    OCPI.settings.tspan = settings.tspan
    OCPI.settings.di = (settings.tspan[2] - settings.tspan[1])/(settings.N-1)

    OCPI.settings.int_method = settings.int_method

    if OCPI.settings.int_method == :Coll
        if !haskey(settings, :Coll_set)
            error("Collocation settings must be provided for Coll integration method.")
        end
        OCPI.settings.Coll_set.order = settings.Coll_set.order
        OCPI.settings.Coll_set.poly = settings.Coll_set.poly
        LGPoly.set_coll_method!(OCPI.settings.Coll_set)
    elseif OCPI.settings.int_method == :EE
        # placeholder
    else
        error("Unsupported integration method: $(OCPI.settings.int_method). Supported methods are :LGR and :EE.")
    end

    OCPI.settings.Params_mode = settings.param_mode

end

function set_x_simple_bounds!(OCPI::OCPInterface, xL::Array{Float64}, xU::Array{Float64})
    @assert length(xL) == OCPI.nx "Length of xL must match number of state variables"
    @assert length(xU) == OCPI.nx "Length of xU must match number of state variables"
    OCPI.xL = xL
    OCPI.xU = xU
end

function set_x_scales!(OCPI::OCPInterface, xScale::Array{Float64})
    @assert length(xScale) == OCPI.nx "Length of xScale must match number of state variables"
    OCPI.xScale = xScale
end

function set_dx_scales!(OCPI::OCPInterface, dxScale::Array{Float64})
    @assert length(dxScale) == OCPI.nx "Length of dxScale must match number of state variables"
    OCPI.dxScale = dxScale
end

function set_u_scales!(OCPI::OCPInterface, uScale::Array{Float64})
    @assert length(uScale) == OCPI.nu "Length of uScale must match number of control inputs"
    OCPI.uScale = uScale
end

function set_p_scales!(OCPI::OCPInterface, pScale::Array{Float64})
    @assert length(pScale) == OCPI.np "Length of pScale must match number of optimization parameters"
    OCPI.pScale = pScale
end

function set_u_simple_bounds!(OCPI::OCPInterface, uL::Array{Float64}, uU::Array{Float64})
    @assert length(uL) == OCPI.nu "Length of uL must match number of control inputs"
    @assert length(uU) == OCPI.nu "Length of uU must match number of control inputs"
    OCPI.uL = uL
    OCPI.uU = uU
end

function set_p_simple_bounds!(OCPI::OCPInterface, pL::Array{Float64}, pU::Array{Float64})
    @assert length(pL) == OCPI.np "Length of pL must match number of optimization parameters"
    @assert length(pU) == OCPI.np "Length of pU must match number of optimization parameters"
    OCPI.pL = pL
    OCPI.pU = pU
end


function load_LUT2d!(OCPI::OCPInterface, filepath::String, LUTName::Symbol; x1_name::String="X", x2_name::String="Y", y_name::String="Z")
    lutData = matread(filepath)
    OCPI.LUTs2D[LUTName] = LUTs.build2D(lutData[x1_name][:], lutData[x2_name][:], lutData[y_name], OCPI, LUTName)
end

function load_LUT1d!(OCPI::OCPInterface, filepath::String, LUTName::Symbol; x_name::String="X", y_name::String="Y")
    lutData = matread(filepath)
    OCPI.LUTs1D[LUTName] = LUTs.build1D(lutData[x_name][:], lutData[y_name][:], OCPI, LUTName)
end



function set_opt_vars!(OCPI::OCPInterface; x0=nothing, u0=nothing, p0=nothing)
    if OCPI.settings.int_method == :EE
        EE_vars!(OCPI; x0=x0, u0=u0, p0=p0)
    elseif OCPI.settings.int_method == :Coll
        Coll_vars!(OCPI; x0=x0, u0=u0, p0=p0)
    else
        error("Unsupported integration method: $(OCPI.settings.int_method). Supported methods are :LGR and :EE.")
    end
end

function EE_vars!(OCPI::OCPInterface; x0=nothing, u0=nothing, p0=nothing)
    if isnothing(x0)
        OCPI.x_n = @variable(OCPI.model, OCPI.xL[i] / OCPI.xScale[i] <= x_n[i in 1:OCPI.nx,1:OCPI.settings.N] <= OCPI.xU[i] / OCPI.xScale[i]) 
    else
        OCPI.x_n = @variable(OCPI.model, OCPI.xL[i] / OCPI.xScale[i] <= x_n[i in 1:OCPI.nx,j in 1:OCPI.settings.N] <= OCPI.xU[i] / OCPI.xScale[i], start = x0[i,j] / OCPI.xScale[i])
        # OCPI.x_n = @variable(OCPI.model, x[i in 1:OCPI.nx,j in 1:OCPI.settings.N], start = x0[i,j] / OCPI.xScale[i])
    end
    OCPI.x = @expression(OCPI.model, [i=1:OCPI.nx, j=1:OCPI.settings.N], OCPI.x_n[i,j] * OCPI.xScale[i])

    if OCPI.nu > 0
        if isnothing(u0)
            OCPI.u_n = @variable(OCPI.model, OCPI.uL[i] / OCPI.uScale[i] <= u_n[i in 1:OCPI.nu,1:OCPI.settings.N-1] <= OCPI.uU[i] / OCPI.uScale[i])
        else
            OCPI.u_n = @variable(OCPI.model, OCPI.uL[i] / OCPI.uScale[i] <= u_n[i in 1:OCPI.nu,j in 1:OCPI.settings.N-1] <= OCPI.uU[i] / OCPI.uScale[i], start = u0[i,j] / OCPI.uScale[i])
            # OCPI.u_n = @variable(OCPI.model, u[i in 1:OCPI.nu,j in 1:OCPI.settings.N-1], start = u0[i,j] / OCPI.uScale[i])
        end
        OCPI.u = @expression(OCPI.model, [i=1:OCPI.nu, j=1:OCPI.settings.N-1], OCPI.u_n[i,j] * OCPI.uScale[i]) 
    end

    if OCPI.np > 0
        if isnothing(p0)
            OCPI.p_n = @variable(OCPI.model, OCPI.pL[i] / OCPI.pScale[i] <= p_n[i=1:OCPI.np] <= OCPI.pU[i] / OCPI.pScale[i]) 
        else
            OCPI.p_n = @variable(OCPI.model, OCPI.pL[i] / OCPI.pScale[i] <= p_n[i=1:OCPI.np] <= OCPI.pU[i] / OCPI.pScale[i], start = p0[i] / OCPI.pScale[i]) 
        end
        OCPI.p = @expression(OCPI.model, [i=1:OCPI.np], OCPI.p_n[i] * OCPI.pScale[i])
    end
end




function Coll_vars!(OCPI::OCPInterface; x0=nothing, u0=nothing, p0=nothing)
    # TODO interpolate initial guess to collocation points
    if isnothing(x0)
        OCPI.x_col_n = @variable(OCPI.model, OCPI.xL[i] / OCPI.xScale[i] <= x_col_n[i in 1:OCPI.nx,j in 1:OCPI.settings.N-1, k in 1:OCPI.settings.Coll_set.order+1] <= OCPI.xU[i] / OCPI.xScale[i]) 
    else
        OCPI.x_col_n = @variable(OCPI.model, OCPI.xL[i] / OCPI.xScale[i] <= x_col_n[i in 1:OCPI.nx,j in 1:OCPI.settings.N-1, k in 1:OCPI.settings.Coll_set.order+1] <= OCPI.xU[i] / OCPI.xScale[i], start = x0[i,j+1] / OCPI.xScale[i])
    end
    OCPI.x_col = @expression(OCPI.model, [i=1:OCPI.nx, j=1:OCPI.settings.N-1, k=1:OCPI.settings.Coll_set.order+1], OCPI.x_col_n[i,j,k] * OCPI.xScale[i]) # j should be 2->N, but this is indexing the vector

    if isnothing(x0)
        x_start_n = @variable(OCPI.model, OCPI.xL[i] / OCPI.xScale[i] <= x_start_n[i in 1:OCPI.nx] <= OCPI.xU[i] / OCPI.xScale[i]) 
    else
        x_start_n = @variable(OCPI.model, OCPI.xL[i] / OCPI.xScale[i] <= x_start_n[i in 1:OCPI.nx] <= OCPI.xU[i] / OCPI.xScale[i], start = x0[i,1] / OCPI.xScale[i])
    end
    x_start = @expression(OCPI.model, [i in 1:OCPI.nx], x_start_n[i] * OCPI.xScale[i])
    x_poly  = @expression(OCPI.model, [i in 1:OCPI.nx, j in 1:OCPI.settings.N-1], sum( OCPI.settings.Coll_set.D .* OCPI.x_col[i,j,:] ) ) # j should be 2->N, but this is indexing the vector

    OCPI.x = hcat(x_start, x_poly)


    if OCPI.nu > 0
        if isnothing(u0)
            OCPI.u_n = @variable(OCPI.model, OCPI.uL[i] / OCPI.uScale[i] <= u_n[i in 1:OCPI.nu,1:OCPI.settings.N-1] <= OCPI.uU[i] / OCPI.uScale[i])
        else
            OCPI.u_n = @variable(OCPI.model, OCPI.uL[i] / OCPI.uScale[i] <= u_n[i in 1:OCPI.nu,j in 1:OCPI.settings.N-1] <= OCPI.uU[i] / OCPI.uScale[i], start = u0[i,j] / OCPI.uScale[i])
        end
        OCPI.u = @expression(OCPI.model, [i=1:OCPI.nu, j=1:OCPI.settings.N-1], OCPI.u_n[i,j] * OCPI.uScale[i]) 
    end

    if OCPI.np > 0
        if OCPI.settings.Params_mode == :single
             if isnothing(p0)
                OCPI.p_n = @variable(OCPI.model, OCPI.pL[i] / OCPI.pScale[i] <= p_n[i=1:OCPI.np] <= OCPI.pU[i] / OCPI.pScale[i]) 
            else
                OCPI.p_n = @variable(OCPI.model, OCPI.pL[i] / OCPI.pScale[i] <= p_n[i=1:OCPI.np] <= OCPI.pU[i] / OCPI.pScale[i], start = p0[i] / OCPI.pScale[i]) 
            end
            OCPI.p = @expression(OCPI.model, [i=1:OCPI.np], OCPI.p_n[i] * OCPI.pScale[i])
        elseif OCPI.settings.Params_mode == :sparse
             if isnothing(p0)
                OCPI.p_n = @variable(OCPI.model, OCPI.pL[i] / OCPI.pScale[i] <= p_n[i=1:OCPI.np, j in 1:OCPI.settings.N-1] <= OCPI.pU[i] / OCPI.pScale[i]) 
            else
                OCPI.p_n = @variable(OCPI.model, OCPI.pL[i] / OCPI.pScale[i] <= p_n[i=1:OCPI.np, j in 1:OCPI.settings.N-1] <= OCPI.pU[i] / OCPI.pScale[i], start = p0[i] / OCPI.pScale[i]) 
            end
            OCPI.p = @expression(OCPI.model, [i=1:OCPI.np, j=1:OCPI.settings.N-1], OCPI.p_n[i,j] * OCPI.pScale[i])
        else
            error("Unsupported parameter mode: $(OCPI.settings.Params_mode). Supported modes are :single and :sparse.")
        end
        
    end
end








function get_param_dict(OCPI::OCPInterface)
    param_dict = copy(ModelingToolkit.initial_conditions(OCPI.sys))
        for var in ModelingToolkit.unknowns(OCPI.sys)
            if haskey(param_dict, var)
                pop!(param_dict, var)
            end
        end
    
    # 1D LUTs
    itp_key = [ k for k in keys(param_dict) 
            if occursin("interpolator", string(k)) && !occursin("2d", string(k)) ]
    for kk in itp_key
        # expr = observed(OCPI.sys)[28].rhs
        name = Symbol(replace(repr(kk.f), "₊interpolator" => ""))

        bad_value = param_dict[kk]
        itp_operator = OCPI.LUTs1D[name].itp_operator
        param_dict = Dict(substitute(k, Dict(kk => kk.f)) => substitute(v, Dict(bad_value => itp_operator)) for (k, v) in param_dict)
    end


    # 2D LUTs
    itp2_key = [ k for k in keys(param_dict) if occursin("interpolator2d", string(k)) ]
    for kk in itp2_key
        # expr = observed(OCPI.sys)[28].rhs
        name = Symbol(replace(repr(kk.f), "₊interpolator2d" => ""))

        bad_value = param_dict[kk]
        itp_operator = OCPI.LUTs2D[name].itp_operator
        param_dict = Dict(substitute(k, Dict(kk => kk.f)) => substitute(v, Dict(bad_value => itp_operator)) for (k, v) in param_dict)
    end

    # write param_dict to txt for debugging
    # open("param_dict.txt", "w") do io
    #     for (k, v) in param_dict
    #         println(io, "$(k) => $(v)")
    #     end
    # end


    
    return param_dict
end


function set_dyn_constraints!(OCPI::OCPInterface;
    f_scale::Union{Float64,Vector{JuMP.AffExpr}, Vector{JuMP.QuadExpr},Vector{JuMP.NonlinearExpr}, Matrix{JuMP.AffExpr}, Matrix{JuMP.QuadExpr}, Matrix{JuMP.NonlinearExpr}}=1.0)
    
    
    if OCPI.settings.int_method == :EE
        EE_dyn_constraints!(OCPI; f_scale=f_scale)
    elseif OCPI.settings.int_method == :Coll
        Coll_dyn_constraints!(OCPI; f_scale=f_scale)
    end
end

function EE_dyn_constraints!(OCPI::OCPInterface;
    f_scale::Union{Float64,Vector{JuMP.AffExpr}, Vector{JuMP.QuadExpr},Vector{JuMP.NonlinearExpr}}=1.0)
        # for now direct transcription with Euler Explicit
        # (x_i+1 - x_i)/dt = f(x_i, u_i, p)
        # (x_i+1 - x_i)/ds = f(x_i, u_i, p) * dt/ds
        # x_i+1 = x_i + f(x_i, u_i, p) * sf * ds
        param_dict = copy(ModelingToolkit.defaults(OCPI.sys))
        for var in ModelingToolkit.unknowns(OCPI.sys)
            if haskey(param_dict, var)
                pop!(param_dict, var)
            end
        end

        # Pre-allocate and compile functions once
        dx = Vector{Any}(undef, OCPI.nx)
        for j in 1:OCPI.nx
            dxj_expr = ModelingToolkit.full_equations(OCPI.sys)[j].rhs
            # Substitute all parameters once, not in a loop
            dxj_expr = SymbolicUtils.substitute(dxj_expr, param_dict)
            dx[j] = build_function(
                dxj_expr,
                OCPI.decision_vars..., 
                expression = Val{false}
            )
        end
        
        OCPI.gDyn.lhs = @expression(OCPI.model, [j in 1:OCPI.nx, i in 1:(OCPI.settings.N-1)], 
            ( OCPI.x[j, i+1] - OCPI.x[j, i] ) / OCPI.dxScale[j] )
        OCPI.gDyn.rhs = @expression(OCPI.model, [j in 1:OCPI.nx, i in 1:(OCPI.settings.N-1)], 
            # ( dx[j](OCPI.x[:, i]..., OCPI.u[:, i]..., OCPI.p_n...) * f_scale[i] ) / OCPI.dxScale[j]) # EE
            ( dx[j](OCPI.x[:, i]..., OCPI.u[:, i]..., OCPI.p...) * f_scale[i] ) / OCPI.dxScale[j]) # EI

        OCPI.gDyn.g = JuMP.@constraint(OCPI.model, [j in 1:OCPI.nx, i in 1:(OCPI.settings.N-1)], 
            OCPI.gDyn.lhs[j,i] == OCPI.gDyn.rhs[j,i])
end

function Coll_dyn_constraints!(OCPI::OCPInterface;
    f_scale::Union{Float64,Matrix{JuMP.AffExpr}, Matrix{JuMP.QuadExpr},Matrix{JuMP.NonlinearExpr}}=1.0)
        param_dict = get_param_dict(OCPI)

        # Pre-allocate and compile functions once
        dx = Vector{Any}(undef, OCPI.nx)
        for j in 1:OCPI.nx
            dxj_expr = ModelingToolkit.full_equations(OCPI.sys)[j].rhs
            # Substitute all parameters once, not in a loop
            dxj_expr = SymbolicUtils.substitute(dxj_expr, param_dict)
            dx[j] = build_function(
                dxj_expr,
                OCPI.decision_vars..., 
                expression = Val{false}
            )
            
        end
        
        OCPI.gDyn.lhs = @expression(OCPI.model, [i in 1:OCPI.nx, j in 1:(OCPI.settings.N-1), k in 2:OCPI.settings.Coll_set.order+1], 
            ( OCPI.settings.Coll_set.C[k,:]' * OCPI.x_col[i,j,:] ) / OCPI.dxScale[i] )

        # rhs = @variable(OCPI.model, [i in 1:OCPI.nx, j in 1:(OCPI.settings.N-1), k in 1:OCPI.settings.Coll_set.order])
        if OCPI.settings.Params_mode == :single || OCPI.np == 0
            OCPI.gDyn.rhs = @expression(OCPI.model, [i in 1:OCPI.nx, j in 1:(OCPI.settings.N-1), k in 2:OCPI.settings.Coll_set.order+1], 
                ( dx[i](OCPI.x_col[:, j, k]..., OCPI.u[:, j]..., OCPI.p...) * f_scale[j,k] ) / OCPI.dxScale[i])
        elseif OCPI.settings.Params_mode == :sparse
            OCPI.gDyn.rhs = @expression(OCPI.model, [i in 1:OCPI.nx, j in 1:(OCPI.settings.N-1), k in 2:OCPI.settings.Coll_set.order+1], 
                ( dx[i](OCPI.x_col[:, j, k]..., OCPI.u[:, j]..., OCPI.p[:, j]...) * f_scale[j,k] ) / OCPI.dxScale[i])
        end
        # JuMP.@constraint(OCPI.model, [i in 1:OCPI.nx, j in 1:(OCPI.settings.N-1), k in 1:OCPI.settings.Coll_set.order], 
        #         rhs[i,j,k] == OCPI.gDyn.rhs[i,j,k])


        OCPI.gDyn.g = JuMP.@constraint(OCPI.model, [i in 1:OCPI.nx, j in 1:(OCPI.settings.N-1), k in 1:OCPI.settings.Coll_set.order], 
            # OCPI.gDyn.lhs[i,j,k] == rhs[i,j,k])
            OCPI.gDyn.lhs[i,j,k] == OCPI.gDyn.rhs[i,j,k])

        # TODO augment state for non LGR collocation to enforce continuity
        OCPI.gDyn.closure = JuMP.@constraint(OCPI.model, [i in 1:OCPI.nx, j in 1:OCPI.settings.N-1], 
            OCPI.x_col[i,j,1] == OCPI.x[i,j] )

end




function set_x0!(OCPI::OCPInterface, x0)
    # @assert size(x0, 1) == OCPI.nx "x0 must have OCPI.nx rows"

    # x0 from dict to vector
    clean_names = [ Symbol(replace(name, "(t)" => "")) for name in OCPI.x_names ]
    x0 = [ x0[Symbol(name)] for name in clean_names ]

    OCPI.gX0 = JuMP.@constraint(OCPI.model, [i in 1:OCPI.nx], 
        OCPI.x[i,1] / OCPI.xScale[i] == x0[i] / OCPI.xScale[i])
end






function set_sparse_param_closure!(OCPI::OCPInterface)
    if OCPI.settings.Params_mode == :sparse
        OCPI.gAux[:param_closure] = JuMP.@constraint(OCPI.model, [i in 1:OCPI.np, j in 1:(OCPI.settings.N-2)], 
            OCPI.p_n[i,j] == OCPI.p_n[i,j+1])
    end
end







function set_y_expr!(OCPI::OCPInterface)
    y_eqs = ModelingToolkit.observed(OCPI.sys)
    param_dict = get_param_dict(OCPI)

    y_exprs = Matrix{Any}(undef, length(y_eqs), OCPI.settings.N-1)
    for i in eachindex(y_eqs)
        y_def = y_eqs[i].rhs
        y_def = SymbolicUtils.substitute(y_def, param_dict)

        y_func = build_function(
            y_def, 
            OCPI.decision_vars..., ModelingToolkit.observables(OCPI.sys)[1:i-1]...,
            expression=Val{false}
        )
        for j in 1:OCPI.settings.N-1
            if OCPI.settings.Params_mode == :single || OCPI.np == 0
                    y_exprs[i, j] = JuMP.@expression(OCPI.model,
                        y_func(OCPI.x[:, j]..., OCPI.u[:, j]..., OCPI.p..., y_exprs[1:i-1, j]...))
            elseif OCPI.settings.Params_mode == :sparse
                    y_exprs[i, j] = JuMP.@expression(OCPI.model,
                        y_func(OCPI.x[:, j]..., OCPI.u[:, j]..., OCPI.p[:, j]..., y_exprs[1:i-1, j]...))
            else
                error("Unsupported parameter mode: $(OCPI.settings.Params_mode). Supported modes are :single and :sparse.")
            end
        end
    end
    OCPI.y_exprs = y_exprs
end



function set_control_acc!(OCPI::OCPInterface; dt::Vector{JuMP.NonlinearExpr})
   N = OCPI.settings.N
    if OCPI.nu > 0
          OCPI.uAcc = @expression(OCPI.model, [i in 1:OCPI.nu, j in 2:(N-2)], 
            (OCPI.u[i,j+1] - 2*OCPI.u[i,j] + OCPI.u[i,j-1]) / (dt[j])^2 )
    end

    delta_idx = findfirst(x -> x == "delta(t)", OCPI.u_names)
    # TODO this is hardcoded for now, need to generalize in settings
    # @constraint(OCPI.model, [i in [delta_idx], j in 1:(N-3)], 
    #     OCPI.uAcc[i,j] <=  1.2)
        # OCPI.uAcc[i,j] <=  3.0)
    # @constraint(OCPI.model, [i in [delta_idx], j in 1:(N-3)], 
    #     OCPI.uAcc[i,j] >= -1.2)
        # OCPI.uAcc[i,j] >= -3.0)
        
end





function get_yi_by_name(OCPI::OCPInterface, name::Symbol)
    y_eqs = ModelingToolkit.observed(OCPI.sys)
    idx = findfirst(eq -> ModelingToolkit.getname(eq.lhs) == name, y_eqs)

    y_expr_vec = OCPI.y_exprs[idx, :]
    return y_expr_vec
end

function collect_x_res!(OCPI::OCPInterface)
    x_res = Dict{Symbol, Vector{Float64}}()
    for i in 1:OCPI.nx
        res_x = value(OCPI.x[i,:])
        name_x = Symbol(replace(OCPI.x_names[i],"(t)"=>""))
        x_res[name_x] = res_x
    end
    OCPI.x_res = x_res
end

function collect_u_res!(OCPI::OCPInterface)
    u_res = Dict{Symbol, Vector{Float64}}()
    for i in 1:OCPI.nu
        res_u = value(OCPI.u[i,:])
        name_u = Symbol(replace(OCPI.u_names[i],"(t)"=>""))
        u_res[name_u] = res_u
    end
    OCPI.u_res = u_res
end

function collect_y_res!(OCPI::OCPInterface)
    y_res = Dict{Symbol, Vector{Float64}}()
    y_eqs = ModelingToolkit.observed(OCPI.sys)
    for i in eachindex(y_eqs)
        name_y = ModelingToolkit.getname(y_eqs[i].lhs)
        res_y = value.(get_yi_by_name(OCPI,name_y))
        y_res[name_y] = res_y
    end
    OCPI.y_res = y_res
end

function collect_p_res!(OCPI::OCPInterface)
    p_res = Dict{Symbol, Vector{Float64}}()
    for i in 1:OCPI.np
        res_p = value(OCPI.p[i,1])
        name_p = Symbol(replace(OCPI.p_names[i],"(t)"=>""))
        p_res[name_p] = [res_p]
    end
    OCPI.p_res = p_res
end


function collect_gDyn_dual!(OCPI::OCPInterface)
    gDyn_dual = Dict{Symbol, Matrix{Float64}}()
    for i in 1:OCPI.nx
        name_x = Symbol(replace(OCPI.x_names[i],"(t)"=>""))
        duals_gDyn = dual.(OCPI.gDyn.g[i,:,:])
        gDyn_dual[name_x] = duals_gDyn
    end
    OCPI.gDyn_dual = gDyn_dual
end

function collect_x_col_res!(OCPI::OCPInterface)
    x_col_res = Dict{Symbol, Matrix{Float64}}()
    for i in 1:OCPI.nx
        res_x_col = value.(OCPI.x_col[i,:,:])
        name_x = Symbol(replace(OCPI.x_names[i],"(t)"=>""))
        x_col_res[name_x] = res_x_col
    end
    OCPI.x_col_res = x_col_res
end

function collect_gAux_dual!(OCPI::OCPInterface)
    gAux_dual = Dict{Symbol, Array{Float64}}()
    for (key, constraints) in OCPI.gAux
        duals_gAux = [ dual(con) for con in constraints ]
        
        if isa(duals_gAux, Vector)
            gAux_dual[key] = reshape(duals_gAux, :, 1)
        else
            gAux_dual[key] = duals_gAux
        end
    end
    OCPI.gAux_dual = gAux_dual
end

function collect_SLS!(OCPI::OCPInterface; f_scale::Union{Float64,Matrix{Float64},Matrix{JuMP.AffExpr},Matrix{JuMP.QuadExpr},Matrix{JuMP.NonlinearExpr}}=1.0)
    N = OCPI.settings.N
    nx = OCPI.nx
    sys = OCPI.sys
    decision_vars = OCPI.decision_vars
    x_names = OCPI.x_names
    u_names = OCPI.u_names
    u_res = OCPI.u_res
    x_col_res = OCPI.x_col_res
    gDyn_dual = OCPI.gDyn_dual
    dxScale = OCPI.dxScale

    gDyn_dual_mat = zeros(nx, N - 1, OCPI.settings.Coll_set.order)
    for i in 1:nx
        name_x = Symbol(replace(x_names[i], "(t)" => ""))
        gDyn_dual_mat[i, :, :] = gDyn_dual[name_x]
    end

    x_col = zeros(nx, N - 1, OCPI.settings.Coll_set.order + 1)
    for i in 1:nx
        name_x = Symbol(replace(x_names[i], "(t)" => ""))
        x_col[i, :, :] = x_col_res[name_x]
    end

    clean_names = [Symbol(replace(string(name), r"\(t\)" => "")) for name in u_names]
    u_col = hcat([u_res[name] for name in clean_names]...)'
    u_col = repeat(u_col, 1, 1, OCPI.settings.Coll_set.order + 1)

    pp = parameters(sys)
    param_dict = get_param_dict(OCPI)

    eqs = full_equations(sys)
    rhs_exprs = [eq.rhs for eq in eqs]

    SLS = Dict{Symbol, Vector{Float64}}()
    for param in pp
        df_dmu_symbolic = Symbolics.derivative.(rhs_exprs, param)
        df_dmu_symbolic = SymbolicUtils.substitute(df_dmu_symbolic, param_dict)

        df_dmu_func = Vector{Any}(undef, nx)
        for i in 1:nx
            df_dmu_func[i] = build_function(df_dmu_symbolic[i], decision_vars...,
                expression=Val{false})
        end


        ls = zeros(N - 1)

        if f_scale == 1.0
            f_scale = ones(N - 1, OCPI.settings.Coll_set.order)
        end


        for i in 1:N-1
            step_sensitivity = 0.0

            for j in 1:OCPI.settings.Coll_set.order
                # 1. Extract states and parameters at the specific collocation point
                x_ij = x_col[:, i, j]
                u_ij = u_col[:, i, j]

                # 2. Evaluate the parameter Jacobian 
                df_dmu_ij = zeros(nx)
                for k in 1:nx
                    df_dmu_ij[k] = df_dmu_func[k](x_ij..., u_ij...)
                end

                # 3. Get the corresponding dual variable vector
                lambda_ij = gDyn_dual_mat[:, i, j] .* dxScale .* f_scale[i, j]


                # 4. Compute the sensitivity for this specific point
                s_ij = dot(lambda_ij, df_dmu_ij)

                # Accumulate
                step_sensitivity += -s_ij
            end

            ls[i] = step_sensitivity
        end

        param_name = Symbol(replace(string(param), r"\(t\)" => ""))
        SLS[Symbol(param_name)] = ls
    end

    OCPI.SLS = SLS
end


function collect_all_results!(OCPI::OCPInterface)
    collect_x_res!(OCPI)
    collect_u_res!(OCPI)
    collect_y_res!(OCPI)
    collect_p_res!(OCPI)
    if OCPI.settings.int_method == :Coll
        collect_x_col_res!(OCPI)
    end
    collect_gDyn_dual!(OCPI)
    collect_gAux_dual!(OCPI)
end


end

