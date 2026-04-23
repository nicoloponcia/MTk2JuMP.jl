
function get_SB()
    
xL = (
    x = -3.0,
    v = -5.0,
    θ = -5*pi,
    ω = -20.0,
    a = -30.0,
)


xU = (
    x =  3.0,
    v =  5.0,
    θ =  5*pi,
    ω =  20.0,
    a =  30.0,
)

uL = (
    F = -10.0,
)
    
uU = (
    F = 10.0,
)

pL = (
    
)


pU = (
    
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
    x = 3.0,
    v = 5.0,
    θ = 5*pi,
    ω = 20.0,
    a = 100.0,
)

u = (
    F = 10.0,
)

dx = (          # state derivative scales
    x = 10.0,
    v = 20.0,
    θ = 2.0,
    ω = 10.0,
    a = 60.0,
)

p = (
    
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