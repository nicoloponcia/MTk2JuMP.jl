@inline _clean_result_name(name::AbstractString) = Symbol(replace(name, "(t)" => ""))
@inline _clean_result_name(name) = Symbol(replace(string(name), "(t)" => ""))

function _collect_named_vector_dict(getter, names)
    result = Dict{Symbol, Vector{Float64}}()
    sizehint!(result, length(names))

    @inbounds for i in eachindex(names)
        result[_clean_result_name(names[i])] = getter(i)
    end

    return result
end

function _collect_named_matrix_dict(getter, names)
    result = Dict{Symbol, Matrix{Float64}}()
    sizehint!(result, length(names))

    @inbounds for i in eachindex(names)
        result[_clean_result_name(names[i])] = getter(i)
    end

    return result
end

function _collect_named_scalar_dict(getter, names)
    result = Dict{Symbol, Vector{Float64}}()
    sizehint!(result, length(names))

    @inbounds for i in eachindex(names)
        result[_clean_result_name(names[i])] = [getter(i)]
    end

    return result
end

function collect_x_res!(OCPI::OCPInterface_)
    OCPI.res.x = _collect_named_vector_dict(i -> value.(view(OCPI.vars.x, i, :)), OCPI.meta.x_names)
end

function collect_u_res!(OCPI::OCPInterface_)
    OCPI.res.u = _collect_named_vector_dict(i -> value.(view(OCPI.vars.u, i, :)), OCPI.meta.u_names)
end

function collect_y_res!(OCPI::OCPInterface_)
    y_eqs = ModelingToolkit.observed(OCPI.sys)
    y_res = Dict{Symbol, Vector{Float64}}()
    sizehint!(y_res, length(y_eqs))

    @inbounds for i in eachindex(y_eqs)
        name_y = _clean_result_name(ModelingToolkit.getname(y_eqs[i].lhs))
        y_res[name_y] = value.(view(OCPI.y_exprs, i, :))
    end

    OCPI.res.y = y_res
end

function collect_p_res!(OCPI::OCPInterface_)
    OCPI.res.p = _collect_named_scalar_dict(i -> value(OCPI.vars.p[i, 1]), OCPI.meta.p_names)
end

function collect_gDyn_dual!(OCPI::OCPInterface_)
    OCPI.duals.gDyn = _collect_named_matrix_dict(i -> dual.(view(OCPI.gDyn.g, i, :, :)), OCPI.meta.x_names)
end

function collect_x_col_res!(OCPI::OCPInterface_)
    OCPI.res.x_col = _collect_named_matrix_dict(i -> value.(view(OCPI.vars.x_col, i, :, :)), OCPI.meta.x_names)
end

function collect_gAux_dual!(OCPI::OCPInterface_)
    gAux_dual = Dict{Symbol, Array{Float64}}()
    sizehint!(gAux_dual, length(OCPI.gAux))

    for (key, constraints) in OCPI.gAux
        duals_gAux = dual.(constraints)
        gAux_dual[key] = duals_gAux isa Vector ? reshape(duals_gAux, :, 1) : duals_gAux
    end

    OCPI.duals.gAux = gAux_dual
end

function collect_SLS!(OCPI::OCPInterface_; f_scale::Union{Float64,Matrix{Float64},Matrix{JuMP.AffExpr},Matrix{JuMP.QuadExpr},Matrix{JuMP.NonlinearExpr}}=1.0)
    N = OCPI.settings.Discretization.N
    nx = OCPI.meta.nx
    sys = OCPI.sys
    decision_vars = OCPI.decision_vars
    x_names = OCPI.meta.x_names
    u_names = OCPI.meta.u_names
    u_res = OCPI.res.u
    x_col_res = OCPI.res.x_col
    gDyn_dual = OCPI.duals.gDyn
    dxScale = OCPI.scales.dx
    order = OCPI.settings.Integration.Collocation.order

    gDyn_dual_mat = zeros(nx, N - 1, order)
    for i in 1:nx
        name_x = Symbol(replace(x_names[i], "(t)" => ""))
        gDyn_dual_mat[i, :, :] = gDyn_dual[name_x]
    end

    x_col = zeros(nx, N - 1, order + 1)
    for i in 1:nx
        name_x = Symbol(replace(x_names[i], "(t)" => ""))
        x_col[i, :, :] = x_col_res[name_x]
    end

    clean_names = [Symbol(replace(string(name), r"\(t\)" => "")) for name in u_names]
    u_col = hcat([u_res[name] for name in clean_names]...)'
    u_col = repeat(u_col, 1, 1, order + 1)

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
        scaled_f = f_scale isa Number ? fill(f_scale, N - 1, order) : f_scale

        for i in 1:N-1
            step_sensitivity = 0.0

            for j in 1:order
                x_ij = x_col[:, i, j]
                u_ij = u_col[:, i, j]

                df_dmu_ij = zeros(nx)
                for k in 1:nx
                    df_dmu_ij[k] = df_dmu_func[k](x_ij..., u_ij...)
                end

                lambda_ij = gDyn_dual_mat[:, i, j] .* dxScale .* scaled_f[i, j]
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
    if OCPI.settings.Integration.int_method == :Coll
        collect_x_col_res!(OCPI)
    end
    collect_gDyn_dual!(OCPI)
    collect_gAux_dual!(OCPI)
end




function collect_solver_info!(OCPI::OCPInterface_; reg::Union{Float64, JuMP.NonlinearExpr}=0.0)
    model = OCPI.model
    solver_info = OCPI.solver_info

    solver_info.termination_status = JuMP.termination_status(model)
    solver_info.solver_time = solve_time(OCPI.model)
    
    solver_info.iterations = barrier_iterations(model)
    solver_info.objective_value = JuMP.result_count(model) > 0 ? JuMP.objective_value(model) : NaN
    solver_info.regularization_value = reg isa Number ? reg : value(reg)
    solver_info.dual_value = dual_objective_value(model)

    return solver_info
end