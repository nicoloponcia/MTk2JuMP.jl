module RocketLandingModel

using ModelingToolkit

function rocket_landing()
    # time
    @independent_variables t
    D = Differential(t)

    # controls (Main, Left, and Right thrusters)
    @discretes F_m(t) F_l(t) F_r(t)

    # Define system parameters
    @parameters begin
        M = 10000.0 # Mass of the rocket (kg)
        I = 50000.0 # Moment of inertia (kg*m^2)
        L = 15.0    # Lever arm of lateral thrusters from COM (m)
        g = 9.81    # Gravity (m/s^2)
    end

    # Define state variables
    @variables begin
        x(t)   # Horizontal position
        y(t)   # Vertical position (altitude)
        v_x(t) # Horizontal velocity
        v_y(t) # Vertical velocity
        θ(t)   # Pitch angle (radians, 0 is pointing straight up)
        ω(t)   # Angular velocity
        a_x(t) # Horizontal acceleration (auxiliary)
        a_y(t) # Vertical acceleration (auxiliary)
        α(t)   # Angular acceleration (auxiliary)
    end

    # Define algebraic/auxiliary variables
    @variables begin
        P(t) # Total thrust proxy (useful for fuel consumption penalty)
        x_base(t)
        y_base(t)
    end

    # Define the equations of motion
    eqs = Equation[
        # Kinematic equations
        D(x)~v_x,
        D(y)~v_y,
        D(θ)~ω,
        D(v_x)~a_x,
        D(v_y)~a_y,
        D(ω)~α,

        # Translational Dynamic equations
        M*a_x~F_m*sin(θ)+(F_l-F_r)*cos(θ),
        M*a_y~F_m*cos(θ)-(F_l-F_r)*sin(θ)-M*g,

        # Rotational Dynamic equations
        I*α~(F_l-F_r)*L,

        # Auxiliary equations
        P~F_m+F_l+F_r,
        x_base~x-L*sin(θ),
        y_base~y-L*cos(θ)
    ]

    # Construct the ODESystem
    sys = ODESystem(eqs, t; name=:rocket_landing)
    return mtkcompile(sys)
end

end