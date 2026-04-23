


function collect_x_res!(OCPI::OCPInterface_)
    x_res = Dict{Symbol, Vector{Float64}}()
    for i in 1:OCPI.meta.nx
        res_x = value(OCPI.vars.x[i,:])
        name_x = Symbol(replace(OCPI.meta.x_names[i],"(t)"=>""))
        x_res[name_x] = res_x
    end
    OCPI.res.x = x_res
end

function collect_u_res!(OCPI::OCPInterface_)
    u_res = Dict{Symbol, Vector{Float64}}()
    for i in 1:OCPI.meta.nu
        res_u = value(OCPI.vars.u[i,:])
        name_u = Symbol(replace(OCPI.meta.u_names[i],"(t)"=>""))
        u_res[name_u] = res_u
    end
    OCPI.res.u = u_res
end

function collect_y_res!(OCPI::OCPInterface_)
    y_res = Dict{Symbol, Vector{Float64}}()
    y_eqs = ModelingToolkit.observed(OCPI.sys)
    for i in eachindex(y_eqs)
        name_y = ModelingToolkit.getname(y_eqs[i].lhs)
        res_y = value.(get_yi_by_name(OCPI,name_y))
        y_res[name_y] = res_y
    end
    OCPI.res.y = y_res
end

function collect_p_res!(OCPI::OCPInterface_)
    p_res = Dict{Symbol, Vector{Float64}}()
    for i in 1:OCPI.meta.np
        res_p = value(OCPI.vars.p[i,1])
        name_p = Symbol(replace(OCPI.meta.p_names[i],"(t)"=>""))
        p_res[name_p] = [res_p]
    end
    OCPI.res.p = p_res
end

function collect_gDyn_dual!(OCPI::OCPInterface_)
    gDyn_dual = Dict{Symbol, Matrix{Float64}}()
    for i in 1:OCPI.meta.nx
        name_x = Symbol(replace(OCPI.meta.x_names[i],"(t)"=>""))
        duals_gDyn = dual.(OCPI.gDyn.g[i,:,:])
        gDyn_dual[name_x] = duals_gDyn
    end
    OCPI.duals.gDyn = gDyn_dual
end

function collect_x_col_res!(OCPI::OCPInterface_)
    x_col_res = Dict{Symbol, Matrix{Float64}}()
    for i in 1:OCPI.meta.nx
        res_x_col = value.(OCPI.vars.x_col[i,:,:])
        name_x = Symbol(replace(OCPI.meta.x_names[i],"(t)"=>""))
        x_col_res[name_x] = res_x_col
    end
    OCPI.res.x_col = x_col_res
end

function collect_gAux_dual!(OCPI::OCPInterface_)
    gAux_dual = Dict{Symbol, Array{Float64}}()
    for (key, constraints) in OCPI.gAux
        duals_gAux = [dual(con) for con in constraints]

        if isa(duals_gAux, Vector)
            gAux_dual[key] = reshape(duals_gAux, :, 1)
        else
            gAux_dual[key] = duals_gAux
        end
    end
    OCPI.duals.gAux = gAux_dual
end

function collect_SLS!(OCPI::OCPInterface_; f_scale::Union{Float64,Matrix{Float64},Matrix{JuMP.AffExpr},Matrix{JuMP.QuadExpr},Matrix{JuMP.NonlinearExpr}}=1.0)
    N = OCPI.settings.N
    nx = OCPI.meta.nx
    sys = OCPI.sys
    decision_vars = OCPI.decision_vars
    x_names = OCPI.meta.x_names
    u_names = OCPI.meta.u_names
    u_res = OCPI.res.u
    x_col_res = OCPI.res.x_col
    gDyn_dual = OCPI.duals.gDyn
    dxScale = OCPI.scales.dx

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
                x_ij = x_col[:, i, j]
                u_ij = u_col[:, i, j]

                df_dmu_ij = zeros(nx)
                for k in 1:nx
                    df_dmu_ij[k] = df_dmu_func[k](x_ij..., u_ij...)
                end

                lambda_ij = gDyn_dual_mat[:, i, j] .* dxScale .* f_scale[i, j]
                s_ij = dot(lambda_ij, df_dmu_ij)
                step_sensitivity += -s_ij
            end

            ls[i] = step_sensitivity
        end

        param_name = Symbol(replace(string(param), r"\(t\)" => ""))
        SLS[Symbol(param_name)] = ls
    end

    OCPI.SLS = SLS
end

function collect_all_results!(OCPI::OCPInterface_)
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




function collect_solver_info!(OCPI::OCPInterface_)
    
    #TODO
    OCPI.solver_info.iterations = 0
    OCPI.solver_info.objective_value = JuMP.objective_value(OCPI.model)
    OCPI.solver_info.regularization_value = 0.0
    OCPI.solver_info.solver_time = 0.0
    OCPI.solver_info.termination_status = JuMP.termination_status(OCPI.model)
    OCPI.solver_info.primals = Dict{Symbol, Vector{Float64}}()
end