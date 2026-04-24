using JSON, JLD2, MAT

function load_LUT2d!(OCPI::OCPInterface_, filepath::String, LUTName::Symbol; x1_name::String="X", x2_name::String="Y", y_name::String="Z")
    lutData = JSON.parsefile(filepath)
    
    # Extract and convert types
    X = Float64.(lutData[x1_name])
    Y = Float64.(lutData[x2_name])
    
    # Convert Vector of Vectors to Matrix
    Z = reduce(hcat, Float64.(row) for row in lutData[y_name])

    OCPI.LUTs2D[LUTName] = build_LUT2d(X, Y, Z, OCPI, LUTName)
end

function load_LUT1d!(OCPI::OCPInterface_, filepath::String, LUTName::Symbol; x_name::String="X", y_name::String="Y")
    lutData = JSON.parsefile(filepath)
    
    X = Float64.(lutData[x_name])
    Y = Float64.(lutData[y_name])
    
    OCPI.LUTs1D[LUTName] = build_LUT1d(X, Y, OCPI, LUTName)
end

# function load_LUT2d!(OCPI, filepath::String, LUTName::Symbol; x1_name::String="X", x2_name::String="Y", y_name::String="Z")
#     # Load the entire JLD2 file as a dictionary
#     lutData = JLD2.load(filepath)
    
#     # JLD2 preserves native Julia types, but we broadcast to Float64 for safety
#     X = Float64.(lutData[x1_name])
#     Y = Float64.(lutData[x2_name])
#     Z = Float64.(lutData[y_name])
    
#     OCPI.LUTs2D[LUTName] = LUTs.build2D(X, Y, Z, OCPI, LUTName)
# end

# function load_LUT1d!(OCPI, filepath::String, LUTName::Symbol; x_name::String="X", y_name::String="Y")
#     lutData = JLD2.load(filepath)
    
#     X = Float64.(lutData[x_name])
#     Y = Float64.(lutData[y_name])
    
#     OCPI.LUTs1D[LUTName] = LUTs.build1D(X, Y, OCPI, LUTName)
# end


# function load_LUT2d!(OCPI, filepath::String, LUTName::Symbol; x1_name::String="X", x2_name::String="Y", y_name::String="Z")
#     lutData = matread(filepath)
    
#     # Use vec() to flatten MATLAB's 1xN or Nx1 matrices into standard Julia vectors
#     X = Float64.(vec(lutData[x1_name]))
#     Y = Float64.(vec(lutData[x2_name]))
    
#     # The Z data remains a 2D matrix
#     Z = Float64.(lutData[y_name])
    
#     OCPI.LUTs2D[LUTName] = LUTs.build2D(X, Y, Z, OCPI, LUTName)
# end

# function load_LUT1d!(OCPI, filepath::String, LUTName::Symbol; x_name::String="X", y_name::String="Y")
#     lutData = matread(filepath)
    
#     X = Float64.(vec(lutData[x_name]))
#     Y = Float64.(vec(lutData[y_name]))
    
#     OCPI.LUTs1D[LUTName] = LUTs.build1D(X, Y, OCPI, LUTName)
# end