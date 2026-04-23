
xL = (
    s1 = -10.0,
    s2 = -20.0,
)


xU = (
    s1 = 10.0,
    s2 = 20.0,
)

uL = (
    u1 = -100.0,
    u2 = 0.0,
)
    
uU = (
    u1 = 100.0,
    u2 = 100.0,
)


pL = (
    p1 = -1.0,
    p2 = -1.0,
)


pU = (
    p1 = 1.0,
    p2 = 1.0,
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
    s1  = 10.0,
    s2  = 20.0,
)

u = (
    u1 = 50.0,
    u2 = 1000.0,
)

dx = (          # state derivative scales
    s1 = 1.0,
    s2 = 1.0,
)

p = (
    p1 = 0.5,
    p2 = 0.5,
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