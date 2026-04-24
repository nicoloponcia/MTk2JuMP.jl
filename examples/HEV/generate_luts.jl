using JSON

function generate_luts()
    # ---------------------------------------------------------
    # 1. Generate Battery 1D LUTs (Voc and R0)
    # ---------------------------------------------------------
    soc_grid = collect(0.0:0.02:1.0) # 0% to 100% SOC

    # Voc: Open Circuit Voltage (V)
    # Realistic Li-ion curve equation
    voc_data = @. 3.2 + 0.8 * soc_grid - 0.2 * soc_grid^2 + 0.1 * exp(-25.0 * soc_grid) - 0.05 * exp(15.0 * (soc_grid - 1.0))

    # R0: Internal Resistance (Ohms)
    # U-shaped resistance curve
    r0_data = @. 0.05 + 0.1 * exp(-20.0 * soc_grid) + 0.08 * exp(-15.0 * (1.0 - soc_grid))

    # Save 1D LUTs to JSON
    voc_dict = Dict("X" => soc_grid, "Y" => voc_data)
    r0_dict = Dict("X" => soc_grid, "Y" => r0_data)

    open("examples/HEV/voc_lut.json", "w") do f
        JSON.print(f, voc_dict, 4)
    end

    open("examples/HEV/r0_lut.json", "w") do f
        JSON.print(f, r0_dict, 4)
    end

    # ---------------------------------------------------------
    # 2. Generate Engine 2D LUT (BSFC)
    # ---------------------------------------------------------
    omega_grid = collect(1000.0:250.0:6000.0) # Engine Speed (RPM)
    torque_grid = collect(10.0:10.0:250.0)    # Engine Torque (Nm)

    # Preallocate matrix
    bsfc_data = zeros(length(omega_grid), length(torque_grid))

    # Define sweet spot
    omega_opt = 2500.0
    torque_opt = 180.0

    for (i, w) in enumerate(omega_grid)
        for (j, t) in enumerate(torque_grid)
            w_norm = (w - omega_opt) / 2000.0
            t_norm = (t - torque_opt) / 100.0
            
            # Quadratic bowl around optimal operating point
            base_bsfc = 220.0 + 40.0 * w_norm^2 + 60.0 * t_norm^2 + 15.0 * w_norm * t_norm
            
            # High penalty for very low torque (pumping losses)
            low_torque_penalty = 150.0 * exp(-0.05 * t)
            
            bsfc_data[i, j] = base_bsfc + low_torque_penalty
        end
    end

    # Save 2D LUT to JSON
    bsfc_dict = Dict(
        "X" => omega_grid,
        "Y" => torque_grid,
        "Z" => bsfc_data
    )

    open("examples/HEV/bsfc_lut.json", "w") do f
        JSON.print(f, bsfc_dict, 4)
    end
end

# Execute the generation
generate_luts()