module CartPoleModel

using ModelingToolkit

function cart_pole()
    # time
    @independent_variables t
    D = Differential(t)

    # controls
    @discretes F(t)

    # Define system parameters
    @parameters begin
        M = 1.0   # Mass of the cart (kg)
        m = 0.1    # Mass of the pole (kg)
        l = 0.5    # Length to the pole's center of mass (m)
        g = 9.81   # Gravity (m/s^2)
    end

    # Define state variables
    @variables begin
        x(t) # Cart position
        v(t) # Cart velocity
        θ(t) # Pole angle (radians, 0 is pointing down)
        ω(t) # Pole angular velocity
        a(t) # Cart acceleration (auxiliary variable)
    end

    # Define algebraic/auxiliary variables (if needed)
    @variables begin
        P(t) # force power
    end

    # Define the equations of motion
    eqs = Equation[
        # Kinematic equations
        D(x)~v,
        D(θ)~ω,
        D(v)~a,

        # Dynamic equations
        (M+m)*a+m*l*D(ω)*cos(θ)-m*l*ω^2*sin(θ)+F~0,
        l*D(ω)+a*cos(θ)+g*sin(θ)~0,

        # power
        P~F*v
    ]

    # Construct the ODESystem
    sys = System(eqs, t; name=:cart_pole)
    return mtkcompile(sys)
end

end