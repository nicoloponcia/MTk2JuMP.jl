module HEVModel

using ModelingToolkit
using MTk2JuMP

function hev(LUTs1D, LUTs2D)
    # time
    @independent_variables t
    D = Differential(t)

    # controls (Engine Torque and Motor Torque)
    @discretes T_e(t) T_m(t)

    # Define system parameters
    @parameters begin
        M_veh = 1500.0   # Vehicle mass (kg)
        r_w = 0.3        # Wheel radius (m)
        FD = 3.5         # Final drive ratio
        Q_bat = 36000.0  # Battery capacity (Amp-seconds or Coulombs)
        rho = 1.2        # Air density (kg/m^3)
        CdA = 0.6        # Drag coefficient * Frontal Area (m^2)
        Crr = 0.01       # Rolling resistance coefficient
        g = 9.81         # Gravity (m/s^2)
    end

    # Define state variables
    @variables begin
        x(t)      # Distance traveled (m)
        v(t)      # Vehicle longitudinal velocity (m/s)
        SOC(t)    # State of Charge (0.0 to 1.0)
        m_fuel(t) # Consumed fuel mass (kg)
    end

    # Define algebraic/auxiliary variables
    @variables begin
        omega_shaft(t) # Rotational speed of the powertrain shaft (rad/s)
        F_trac(t)      # Tractive force at the wheels (N)
        F_drag(t)      # Aerodynamic drag force (N)
        F_roll(t)      # Rolling resistance force (N)
        P_batt(t)      # Battery electrical power demand (W)
        I_batt(t)      # Battery current (A)
        V_oc(t)        # Open-circuit voltage (V)
        R_0(t)         # Internal resistance (Ohms)
        mdot_fuel(t)   # Fuel mass flow rate (kg/s)
    end

    @named bsfc_lut = MTk2JuMP.MTk.Interpolations.Interpolation2d(; itp=LUTs2D[:bsfc_lut].itp)
    @named voc_lut = MTk2JuMP.MTk.Interpolations.Interpolation1d(; itp=LUTs1D[:voc_lut].itp)
    @named r0_lut = MTk2JuMP.MTk.Interpolations.Interpolation1d(; itp=LUTs1D[:r0_lut].itp)

    # Define the equations of motion and system dynamics
    eqs = Equation[
        # --- Kinematic & Longitudinal Dynamics ---
        D(x) ~ v,
        F_drag ~ 0.5 * rho * CdA * v^2,
        F_roll ~ M_veh * g * Crr,
        
        # Total tractive force mapped from shaft torques to the wheels
        F_trac ~ ((T_e + T_m) * FD) / r_w,
        D(v) ~ (F_trac - F_drag - F_roll) / M_veh,

        # Powertrain kinematics (assuming wheels don't slip)
        omega_shaft ~ (v / r_w) * FD,

        # --- 2D LUT: Engine Fuel Consumption ---
        bsfc_lut.input1.u ~ omega_shaft,
        bsfc_lut.input2.u ~ T_e,
        mdot_fuel ~ bsfc_lut.output.u,
        # mdot_fuel ~ 1.0,
        D(m_fuel) ~ mdot_fuel,

        # --- 1D LUTs: Battery Characteristics ---
        voc_lut.input.u ~ SOC,
        V_oc ~ voc_lut.output.u,
        r0_lut.input.u ~ SOC,
        R_0 ~ r0_lut.output.u,

        # --- Battery Equivalent Circuit Model ---
        # Assuming 100% motor efficiency for simplicity (you could add another 2D LUT here)
        P_batt ~ T_m * omega_shaft, 
        
        # Solving the circuit equation: P = V_oc * I - I^2 * R_0 using the quadratic formula
        I_batt ~ (V_oc - sqrt(V_oc^2 - 4 * R_0 * P_batt)) / (2 * R_0),
        
        # State of charge depletion
        D(SOC) ~ -I_batt / Q_bat
    ]

    # Construct the ODESystem
    sys = System(eqs, t; name=:hev, systems=[bsfc_lut, voc_lut, r0_lut])
    return mtkcompile(sys)
end

end