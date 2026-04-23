function get_yi_by_name(OCPI::OCPInterface_, name::Symbol)
    y_eqs = ModelingToolkit.observed(OCPI.sys)
    idx = findfirst(eq -> ModelingToolkit.getname(eq.lhs) == name, y_eqs)

    y_expr_vec = OCPI.y_exprs[idx, :]
    return y_expr_vec
end


function set_control_der!(OCPI::OCPInterface_; dt::Union{Vector{JuMP.NonlinearExpr}, Vector{Float64}})
    N = OCPI.settings.N
    dt = length(dt) == 1 && typeof(dt[1]) == Float64 ? fill(dt[1], OCPI.settings.N-1) : dt

    if OCPI.meta.nu > 0
        OCPI.uDer = @expression(OCPI.model, [i in 1:OCPI.meta.nu, j in 1:(N-2)],
            (OCPI.vars.u[i,j+1] - OCPI.vars.u[i,j] ) / (dt[j]))
    end
end



function set_control_acc!(OCPI::OCPInterface_; dt::Union{Vector{JuMP.NonlinearExpr}, Vector{Float64}})
    N = OCPI.settings.N
    dt = length(dt) == 1 && typeof(dt[1]) == Float64 ? fill(dt[1], OCPI.settings.N-2) : dt
        
    if OCPI.meta.nu > 0
        OCPI.uAcc = @expression(OCPI.model, [i in 1:OCPI.meta.nu, j in 2:(N-2)],
            (OCPI.vars.u[i,j+1] - 2*OCPI.vars.u[i,j] + OCPI.vars.u[i,j-1]) / (dt[j])^2)
    end
end