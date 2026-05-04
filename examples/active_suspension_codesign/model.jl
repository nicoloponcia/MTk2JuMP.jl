module ActiveSuspensionModel

using ModelingToolkit


function active_suspension()
    # time
    @independent_variables t
    D = Differential(t)

    # controls and external disturbances
    @discretes F(t) z_r(t) z_r_dot(t)


    # Define system parameters
    @parameters begin
        m = 200.0 # Vehicle quarter-mass (kg)
        k      # Spring stiffness (N/m) [Optimizable parameter]
        c      # Damping coefficient (N*s/m) [Optimizable parameter]
        g = 9.81 # Gravitational acceleration (m/s^2)
    end

    # Define state variables
    @variables begin
        z(t) # Mass displacement
        z_dot(t) # Mass velocity
        a(t)   # Vertical acceleration (auxiliary variable for cost function)
    end

    # Define the equations of motion
    eqs = Equation[
        # Kinematic equations
        D(z) ~ z_dot,
        D(z_dot) ~ a,
        
        a ~ (F - k * (z - z_r) - c * (z_dot - z_r_dot)) / m - g,
    ]

    # Construct the System
    sys = System(eqs, t; name=:active_suspension)
    return mtkcompile(sys)
end

end