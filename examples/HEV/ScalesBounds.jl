
function get_SB()
    
xL = (
    SOC = 0.0,
    m_fuel = -1.0,
    v = 0.0,
    x = 0.0,
)


xU = (
    SOC = 1.0,
    m_fuel = 100.0,
    v = 50.0,
    x = 100000.0,
)

uL = (
    T_e = 0.0,
    T_m = -50.0,
    F_brk = 0.0,
)
    
uU = (
    T_e = 250.0,
    T_m = 50.0,
    F_brk = 10000.0,
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
    SOC = 1.0,
    m_fuel = 100.0,
    v = 50.0,
    x = 100000.0,
)

u = (
    T_e = 250.0,
    T_m = 50.0,
    F_brk = 10000.0,
)

dx = (          # state derivative scales
    SOC = 1.0,
    m_fuel = 1.0,
    v = 20.0,
    x = 50.0,
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