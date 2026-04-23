
function get_SB()
    
xL = (
    ω = -20.0,
    v_y = -50.0,
    v_x = -50.0,
    θ = -pi/2,
    y = 0.0,
    x = -200.0,
)


xU = (
    ω = 20.0,
    v_y = 50.0,
    v_x = 50.0,
    θ = pi/2,
    y = 150.0,
    x = 200.0,
)

uL = (
    F_m = 0.0,
    F_l = -10000.0,
    F_r = -10000.0,
)
    
uU = (
    F_m = 300000.0,
    F_l = 10000.0,
    F_r = 10000.0,
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
    ω = 20.0,
    v_y = 50.0,
    v_x = 50.0,
    θ = pi/2,
    y = 150.0,
    x = 200.0,
)

u = (
    F_m = 100000.0,
    F_l = 10000.0,
    F_r = 10000.0,
)

dx = (          # state derivative scales
    ω = 20.0,
    v_y = 5.0,
    v_x = 5.0,
    θ = 5*pi,
    y = 3.0,
    x = 3.0,
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