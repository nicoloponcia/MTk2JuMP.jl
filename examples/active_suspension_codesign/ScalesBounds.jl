
function get_SB()
    
xL = (
    z = -50e-3,
    z_dot = -10.0,
)


xU = (
    z =  50e-3,
    z_dot = 10.0,
)

# put also input knowns
uL = (
    F = -1000.0,
    z_r = -1, # irrelevant since fixed
    z_r_dot = -1, # irrelevant since fixed
)
    
uU = (
    F = 1000.0,
    z_r = 1, # irrelevant since fixed
    z_r_dot = 1, # irrelevant since fixed
)

pL = (
    k = 0.0,
    c = 0.0,
)


pU = (
    k = 40000.0,
    c = 5000.0,
)


Bounds = (
    xL = xL,
    xU = xU,
    uL = uL,
    uU = uU,
    pL = pL,
    pU = pU,
)


x = (
    z = 1.0,
    z_dot = 1.0,
)

# put also input knowns
u = (
    F = 1000.0,
    z_r = 1.0, # irrelevant since fixed
    z_r_dot = 1.0, # irrelevant since fixed
)

dx = (          # state derivative scales
    z = 1.0,
    z_dot = 10.0,
)

p = (
    k = 40000.0,
    c = 5000.0,
)


Scales = (
    x = x,
    u = u,
    dx = dx,
    p = p,
)


SB = (
    Bounds = Bounds,
    Scales = Scales,
)

return SB
end

SB = get_SB()