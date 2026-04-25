
function get_SB()
    
xL = (
    z = 0.0,
    z_dot = -10.0,
)


xU = (
    z =  1.0,
    z_dot = 10.0,
)

uL = (
    F = -10.0,
)
    
uU = (
    F = 10.0,
)

pL = (
    k = 1000.0,
    c = 100.0,
)


pU = (
    k = 5000.0,
    c = 500.0,
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
    z = 0.5,
    z_dot = 0.0,
)

u = (
    F = 10.0,
)

dx = (          # state derivative scales
    z = 1.0,
    z_dot = 10.0,
)

p = (
    k = 1000.0,
    c = 100.0,
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