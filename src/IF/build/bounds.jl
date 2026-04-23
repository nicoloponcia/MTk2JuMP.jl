function set_bounds!(OCPI::OCPInterface_, Bounds)
    # vectorize bound
    xL, xU = compute_bounds(OCPI.meta.x_names, Bounds.xL, Bounds.xU)
    uL, uU = compute_bounds(OCPI.meta.u_names, Bounds.uL, Bounds.uU)
    pL, pU = compute_bounds(OCPI.meta.p_names, Bounds.pL, Bounds.pU)

    # set bound
    set_x_simple_bounds!(OCPI, xL, xU)
    set_u_simple_bounds!(OCPI, uL, uU)
    set_p_simple_bounds!(OCPI, pL, pU)
end


function compute_bounds(names, L, U)
    # x_names = OCPI.x_names

    lb = Array{Float64}(undef, length(names))
    ub = Array{Float64}(undef, length(names))

    for (i, name) in enumerate(names)
        # Extract base name (e.g., "omega_FL(t)" -> :omega_FL)
        clean_symbol = Symbol(split(name, '(')[1])

        try
            lb[i] = getfield(L, clean_symbol)
            ub[i] = getfield(U, clean_symbol)
        catch
            error("Unknown state variable: $name")
        end
    end

    return lb, ub
end

function set_x_simple_bounds!(OCPI::OCPInterface_, xL::Array{Float64}, xU::Array{Float64})
    @assert length(xL) == OCPI.meta.nx "Length of xL must match number of state variables"
    @assert length(xU) == OCPI.meta.nx "Length of xU must match number of state variables"
    OCPI.bounds.xL = xL
    OCPI.bounds.xU = xU
end

function set_u_simple_bounds!(OCPI::OCPInterface_, uL::Array{Float64}, uU::Array{Float64})
    @assert length(uL) == OCPI.meta.nu "Length of uL must match number of control inputs"
    @assert length(uU) == OCPI.meta.nu "Length of uU must match number of control inputs"
    OCPI.bounds.uL = uL
    OCPI.bounds.uU = uU
end

function set_p_simple_bounds!(OCPI::OCPInterface_, pL::Array{Float64}, pU::Array{Float64})
    @assert length(pL) == OCPI.meta.np "Length of pL must match number of optimization parameters"
    @assert length(pU) == OCPI.meta.np "Length of pU must match number of optimization parameters"
    OCPI.bounds.pL = pL
    OCPI.bounds.pU = pU
end



function set_scales!(OCPI::OCPInterface_, Scales)
    # vectorize scale
    xScale = compute_scales(OCPI.meta.x_names, Scales.x)
    uScale = compute_scales(OCPI.meta.u_names, Scales.u)
    pScale = compute_scales(OCPI.meta.p_names, Scales.p)
    dxScale = compute_scales(OCPI.meta.x_names, Scales.dx)

    # set scale
    set_x_scales!(OCPI, xScale)
    set_u_scales!(OCPI, uScale)
    set_p_scales!(OCPI, pScale)
    set_dx_scales!(OCPI, dxScale)
end


function compute_scales(names, Scales)
    scale_vec = Array{Float64}(undef, length(names))
    for (i, name) in enumerate(names)
        # Extract base name (e.g., "omega_FL(t)" -> :omega_FL)
        clean_symbol = Symbol(split(name, '(')[1])

        try
            scale_vec[i] = getfield(Scales, clean_symbol)
        catch
            error("Unknown variable for scaling: $name")
        end
    end
    return scale_vec
end



function set_x_scales!(OCPI::OCPInterface_, xScale::Array{Float64})
    @assert length(xScale) == OCPI.meta.nx "Length of xScale must match number of state variables"
    OCPI.scales.x = xScale
end

function set_dx_scales!(OCPI::OCPInterface_, dxScale::Array{Float64})
    @assert length(dxScale) == OCPI.meta.nx "Length of dxScale must match number of state variables"
    OCPI.scales.dx = dxScale
end

function set_u_scales!(OCPI::OCPInterface_, uScale::Array{Float64})
    @assert length(uScale) == OCPI.meta.nu "Length of uScale must match number of control inputs"
    OCPI.scales.u = uScale
end

function set_p_scales!(OCPI::OCPInterface_, pScale::Array{Float64})
    @assert length(pScale) == OCPI.meta.np "Length of pScale must match number of optimization parameters"
    OCPI.scales.p = pScale
end