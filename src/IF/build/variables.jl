function set_opt_vars!(OCPI::OCPInterface_; x0=nothing, u0=nothing, p0=nothing)
    if OCPI.settings.int_method == :EE
        EE_vars!(OCPI; x0=x0, u0=u0, p0=p0)
    elseif OCPI.settings.int_method == :Coll
        Coll_vars!(OCPI; x0=x0, u0=u0, p0=p0)
    else
        error("Unsupported integration method: $(OCPI.settings.int_method). Supported methods are :Coll and :EE.")
    end
end

function EE_vars!(OCPI::OCPInterface_; x0=nothing, u0=nothing, p0=nothing)
    if isnothing(x0)
        OCPI.vars_n.x = @variable(OCPI.model, OCPI.bounds.xL[i] / OCPI.scales.x[i] <= x_n[i in 1:OCPI.meta.nx,1:OCPI.settings.N] <= OCPI.bounds.xU[i] / OCPI.scales.x[i])
    else
        OCPI.vars_n.x = @variable(OCPI.model, OCPI.bounds.xL[i] / OCPI.scales.x[i] <= x_n[i in 1:OCPI.meta.nx,j in 1:OCPI.settings.N] <= OCPI.bounds.xU[i] / OCPI.scales.x[i], start = x0[i,j] / OCPI.scales.x[i])
    end
    OCPI.vars.x = @expression(OCPI.model, [i=1:OCPI.meta.nx, j=1:OCPI.settings.N], OCPI.vars_n.x[i,j] * OCPI.scales.x[i])

    if OCPI.meta.nu > 0
        if isnothing(u0)
            OCPI.vars_n.u = @variable(OCPI.model, OCPI.bounds.uL[i] / OCPI.scales.u[i] <= u_n[i in 1:OCPI.meta.nu,1:OCPI.settings.N-1] <= OCPI.bounds.uU[i] / OCPI.scales.u[i])
        else
            OCPI.vars_n.u = @variable(OCPI.model, OCPI.bounds.uL[i] / OCPI.scales.u[i] <= u_n[i in 1:OCPI.meta.nu,j in 1:OCPI.settings.N-1] <= OCPI.bounds.uU[i] / OCPI.scales.u[i], start = u0[i,j] / OCPI.scales.u[i])
        end
        OCPI.vars.u = @expression(OCPI.model, [i=1:OCPI.meta.nu, j=1:OCPI.settings.N-1], OCPI.vars_n.u[i,j] * OCPI.scales.u[i])
    end

    if OCPI.meta.np > 0
        if isnothing(p0)
            OCPI.vars_n.p = @variable(OCPI.model, OCPI.bounds.pL[i] / OCPI.scales.p[i] <= p_n[i=1:OCPI.meta.np] <= OCPI.bounds.pU[i] / OCPI.scales.p[i])
        else
            OCPI.vars_n.p = @variable(OCPI.model, OCPI.bounds.pL[i] / OCPI.scales.p[i] <= p_n[i=1:OCPI.meta.np] <= OCPI.bounds.pU[i] / OCPI.scales.p[i], start = p0[i] / OCPI.scales.p[i])
        end
        OCPI.vars.p = @expression(OCPI.model, [i=1:OCPI.meta.np], OCPI.vars_n.p[i] * OCPI.scales.p[i])
    end
end

function Coll_vars!(OCPI::OCPInterface_; x0=nothing, u0=nothing, p0=nothing)
    if isnothing(x0)
        OCPI.vars_n.x_col = @variable(OCPI.model, OCPI.bounds.xL[i] / OCPI.scales.x[i] <= x_col_n[i in 1:OCPI.meta.nx,j in 1:OCPI.settings.N-1, k in 1:OCPI.settings.Coll_set.order+1] <= OCPI.bounds.xU[i] / OCPI.scales.x[i])
    else
        OCPI.vars_n.x_col = @variable(OCPI.model, OCPI.bounds.xL[i] / OCPI.scales.x[i] <= x_col_n[i in 1:OCPI.meta.nx,j in 1:OCPI.settings.N-1, k in 1:OCPI.settings.Coll_set.order+1] <= OCPI.bounds.xU[i] / OCPI.scales.x[i], start = x0[i,j+1] / OCPI.scales.x[i])
    end
    OCPI.vars.x_col = @expression(OCPI.model, [i=1:OCPI.meta.nx, j=1:OCPI.settings.N-1, k=1:OCPI.settings.Coll_set.order+1], OCPI.vars_n.x_col[i,j,k] * OCPI.scales.x[i])

    if isnothing(x0)
        x_start_n = @variable(OCPI.model, OCPI.bounds.xL[i] / OCPI.scales.x[i] <= x_start_n[i in 1:OCPI.meta.nx] <= OCPI.bounds.xU[i] / OCPI.scales.x[i])
    else
        x_start_n = @variable(OCPI.model, OCPI.bounds.xL[i] / OCPI.scales.x[i] <= x_start_n[i in 1:OCPI.meta.nx] <= OCPI.bounds.xU[i] / OCPI.scales.x[i], start = x0[i,1] / OCPI.scales.x[i])
    end
    x_start = @expression(OCPI.model, [i in 1:OCPI.meta.nx], x_start_n[i] * OCPI.scales.x[i])
    x_poly  = @expression(OCPI.model, [i in 1:OCPI.meta.nx, j in 1:OCPI.settings.N-1], sum(OCPI.settings.Coll_set.D .* OCPI.vars.x_col[i,j,:]))

    OCPI.vars.x = hcat(x_start, x_poly)

    if OCPI.meta.nu > 0
        if isnothing(u0)
            OCPI.vars_n.u = @variable(OCPI.model, OCPI.bounds.uL[i] / OCPI.scales.u[i] <= u_n[i in 1:OCPI.meta.nu,1:OCPI.settings.N-1] <= OCPI.bounds.uU[i] / OCPI.scales.u[i])
        else
            OCPI.vars_n.u = @variable(OCPI.model, OCPI.bounds.uL[i] / OCPI.scales.u[i] <= u_n[i in 1:OCPI.meta.nu,j in 1:OCPI.settings.N-1] <= OCPI.bounds.uU[i] / OCPI.scales.u[i], start = u0[i,j] / OCPI.scales.u[i])
        end
        OCPI.vars.u = @expression(OCPI.model, [i=1:OCPI.meta.nu, j=1:OCPI.settings.N-1], OCPI.vars_n.u[i,j] * OCPI.scales.u[i])
    end

    if OCPI.meta.np > 0
        if OCPI.settings.Params_mode == :single
            if isnothing(p0)
                OCPI.vars_n.p = @variable(OCPI.model, OCPI.bounds.pL[i] / OCPI.scales.p[i] <= p_n[i=1:OCPI.meta.np] <= OCPI.bounds.pU[i] / OCPI.scales.p[i])
            else
                OCPI.vars_n.p = @variable(OCPI.model, OCPI.bounds.pL[i] / OCPI.scales.p[i] <= p_n[i=1:OCPI.meta.np] <= OCPI.bounds.pU[i] / OCPI.scales.p[i], start = p0[i] / OCPI.scales.p[i])
            end
            OCPI.vars.p = @expression(OCPI.model, [i=1:OCPI.meta.np], OCPI.vars_n.p[i] * OCPI.scales.p[i])
        elseif OCPI.settings.Params_mode == :sparse
            if isnothing(p0)
                OCPI.vars_n.p = @variable(OCPI.model, OCPI.bounds.pL[i] / OCPI.scales.p[i] <= p_n[i=1:OCPI.meta.np, j in 1:OCPI.settings.N-1] <= OCPI.bounds.pU[i] / OCPI.scales.p[i])
            else
                OCPI.vars_n.p = @variable(OCPI.model, OCPI.bounds.pL[i] / OCPI.scales.p[i] <= p_n[i=1:OCPI.meta.np, j in 1:OCPI.settings.N-1] <= OCPI.bounds.pU[i] / OCPI.scales.p[i], start = p0[i] / OCPI.scales.p[i])
            end
            OCPI.vars.p = @expression(OCPI.model, [i=1:OCPI.meta.np, j=1:OCPI.settings.N-1], OCPI.vars_n.p[i,j] * OCPI.scales.p[i])
        else
            error("Unsupported parameter mode: $(OCPI.settings.Params_mode). Supported modes are :single and :sparse.")
        end
    end
end
