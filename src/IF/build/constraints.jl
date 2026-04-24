function get_param_dict(OCPI::OCPInterface_)
    param_dict = copy(ModelingToolkit.initial_conditions(OCPI.sys))
    for var in ModelingToolkit.unknowns(OCPI.sys)
        if haskey(param_dict, var)
            pop!(param_dict, var)
        end
    end

    itp_key = [k for k in keys(param_dict)
               if occursin("interpolator", string(k)) && !occursin("2d", string(k))]
    for kk in itp_key
        name = replace(repr(kk.name), "₊interpolator1d" => "")
        name = replace(name, ":" => "")
        name = Symbol(name)

        bad_value = param_dict[kk]
        itp_operator = OCPI.LUTs1D[name].itp_operator
        param_dict = Dict(substitute(k, Dict(kk => kk)) => substitute(v, Dict(bad_value => itp_operator)) for (k, v) in param_dict)
    end

    itp2_key = [k for k in keys(param_dict) if occursin("interpolator2d", string(k))]
    for kk in itp2_key
        name = replace(repr(kk.name), "₊interpolator2d" => "")
        name = replace(name, ":" => "")
        name = Symbol(name)

        bad_value = param_dict[kk]
        itp_operator = OCPI.LUTs2D[name].itp_operator
        param_dict = Dict(substitute(k, Dict(kk => kk)) => substitute(v, Dict(bad_value => itp_operator)) for (k, v) in param_dict)
    end

    return param_dict
end

function set_gDyn!(OCPI::OCPInterface_;
    f_scale::Union{Float64,Vector{Float64}, Matrix{Float64}, JuMP.AbstractJuMPScalar, Vector{<:JuMP.AbstractJuMPScalar}, Matrix{<:JuMP.AbstractJuMPScalar}}=1.0)

    if OCPI.settings.int_method == :EE
        f_scale = (f_scale isa Union{Number, JuMP.AbstractJuMPScalar}) ? fill(f_scale, OCPI.settings.N-1) : f_scale
        EE_dyn_constraints!(OCPI; f_scale=f_scale)
    elseif OCPI.settings.int_method == :Coll
        f_scale = (f_scale isa Union{Number, JuMP.AbstractJuMPScalar}) ? fill(f_scale, OCPI.settings.N-1, OCPI.settings.Coll_set.order+1) : f_scale
        Coll_dyn_constraints!(OCPI; f_scale=f_scale)
    end
end

function EE_dyn_constraints!(OCPI::OCPInterface_;
    f_scale::Union{Float64,Vector{Float64}, Matrix{Float64}, JuMP.AbstractJuMPScalar, Vector{<:JuMP.AbstractJuMPScalar}, Matrix{<:JuMP.AbstractJuMPScalar}}=1.0)
    param_dict = get_param_dict(OCPI)

    dx = Vector{Any}(undef, OCPI.meta.nx)
    for j in 1:OCPI.meta.nx
        dxj_expr = ModelingToolkit.full_equations(OCPI.sys)[j].rhs
        dxj_expr = SymbolicUtils.substitute(dxj_expr, param_dict)
        dx[j] = build_function(
            dxj_expr,
            OCPI.decision_vars...,
            expression = Val{false}
        )
    end

    OCPI.gDyn.lhs = @expression(OCPI.model, [j in 1:OCPI.meta.nx, i in 1:(OCPI.settings.N-1)],
        (OCPI.vars.x[j, i+1] - OCPI.vars.x[j, i]) / OCPI.scales.dx[j])
    OCPI.gDyn.rhs = @expression(OCPI.model, [j in 1:OCPI.meta.nx, i in 1:(OCPI.settings.N-1)],
        (dx[j](OCPI.vars.x[:, i]..., OCPI.vars.u[:, i]..., OCPI.vars.p...) * f_scale[i]) / OCPI.scales.dx[j])

    OCPI.gDyn.g = JuMP.@constraint(OCPI.model, [j in 1:OCPI.meta.nx, i in 1:(OCPI.settings.N-1)],
        OCPI.gDyn.lhs[j,i] == OCPI.gDyn.rhs[j,i])
end

function Coll_dyn_constraints!(OCPI::OCPInterface_;
    f_scale::Union{Float64,Vector{Float64}, Matrix{Float64}, JuMP.AbstractJuMPScalar, Vector{<:JuMP.AbstractJuMPScalar}, Matrix{<:JuMP.AbstractJuMPScalar}}=1.0)
    param_dict = get_param_dict(OCPI)

    dx = Vector{Any}(undef, OCPI.meta.nx)
    for j in 1:OCPI.meta.nx
        dxj_expr = ModelingToolkit.full_equations(OCPI.sys)[j].rhs
        dxj_expr = SymbolicUtils.substitute(dxj_expr, param_dict)
        dx[j] = build_function(
            dxj_expr,
            OCPI.decision_vars...,
            expression = Val{false}
        )
    end

    OCPI.gDyn.lhs = @expression(OCPI.model, [i in 1:OCPI.meta.nx, j in 1:(OCPI.settings.N-1), k in 2:OCPI.settings.Coll_set.order+1],
        (OCPI.settings.Coll_set.C[k,:]' * OCPI.vars.x_col[i,j,:]) / OCPI.scales.dx[i])

    if OCPI.settings.Params_mode == :single || OCPI.meta.np == 0
        OCPI.gDyn.rhs = @expression(OCPI.model, [i in 1:OCPI.meta.nx, j in 1:(OCPI.settings.N-1), k in 2:OCPI.settings.Coll_set.order+1],
            (dx[i](OCPI.vars.x_col[:, j, k]..., OCPI.vars.u[:, j]..., OCPI.vars.p...) * f_scale[j,k]) / OCPI.scales.dx[i])
    elseif OCPI.settings.Params_mode == :sparse
        OCPI.gDyn.rhs = @expression(OCPI.model, [i in 1:OCPI.meta.nx, j in 1:(OCPI.settings.N-1), k in 2:OCPI.settings.Coll_set.order+1],
            (dx[i](OCPI.vars.x_col[:, j, k]..., OCPI.vars.u[:, j]..., OCPI.vars.p[:, j]...) * f_scale[j,k]) / OCPI.scales.dx[i])
    end

    OCPI.gDyn.g = JuMP.@constraint(OCPI.model, [i in 1:OCPI.meta.nx, j in 1:(OCPI.settings.N-1), k in 1:OCPI.settings.Coll_set.order],
        OCPI.gDyn.lhs[i,j,k] == OCPI.gDyn.rhs[i,j,k])

    start_node_idx = findfirst(tau -> isapprox(tau, 0.0; atol=1e-12), OCPI.settings.Coll_set.nodes)
    isnothing(start_node_idx) && error("Collocation nodes must include tau=0.0 to enforce interval closure.")

    closure_init = JuMP.@constraint(OCPI.model, [i in 1:OCPI.meta.nx],
        OCPI.vars.x_col[i,1,start_node_idx] == OCPI.vars.x[i,1])

    if OCPI.settings.N > 2
        closure_link = JuMP.@constraint(OCPI.model, [i in 1:OCPI.meta.nx, j in 2:(OCPI.settings.N-1)],
            OCPI.vars.x_col[i,j,start_node_idx] == sum(OCPI.settings.Coll_set.D .* OCPI.vars.x_col[i,j-1,:]))
        OCPI.gDyn.closure = vcat(vec(closure_init), vec(closure_link))
    else
        OCPI.gDyn.closure = vec(closure_init)
    end

end

function set_x0!(OCPI::OCPInterface_, x0)
    clean_names = [Symbol(replace(name, "(t)" => "")) for name in OCPI.meta.x_names]
    x0 = [x0[Symbol(name)] for name in clean_names]

    OCPI.gX0 = JuMP.@constraint(OCPI.model, [i in 1:OCPI.meta.nx],
        OCPI.vars.x[i,1] / OCPI.scales.x[i] == x0[i] / OCPI.scales.x[i])
end

function set_sparse_param_closure!(OCPI::OCPInterface_)
    if OCPI.settings.Params_mode == :sparse
        OCPI.gAux[:param_closure] = JuMP.@constraint(OCPI.model, [i in 1:OCPI.meta.np, j in 1:(OCPI.settings.N-2)],
            OCPI.vars_n.p[i,j] == OCPI.vars_n.p[i,j+1])
    end
end

function set_y_expr!(OCPI::OCPInterface_)
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
            if OCPI.settings.Params_mode == :single || OCPI.meta.np == 0
                y_exprs[i, j] = JuMP.@expression(OCPI.model,
                    y_func(OCPI.vars.x[:, j]..., OCPI.vars.u[:, j]..., OCPI.vars.p..., y_exprs[1:i-1, j]...))
            elseif OCPI.settings.Params_mode == :sparse
                y_exprs[i, j] = JuMP.@expression(OCPI.model,
                    y_func(OCPI.vars.x[:, j]..., OCPI.vars.u[:, j]..., OCPI.vars.p[:, j]..., y_exprs[1:i-1, j]...))
            else
                error("Unsupported parameter mode: $(OCPI.settings.Params_mode). Supported modes are :single and :sparse.")
            end
        end
    end
    OCPI.y_exprs = y_exprs
end


