function load_LUT2d!(OCPI::OCPInterface_, filepath::String, LUTName::Symbol; x1_name::String="X", x2_name::String="Y", y_name::String="Z")
    lutData = matread(filepath)
    OCPI.LUTs2D[LUTName] = LUTs.build2D(lutData[x1_name][:], lutData[x2_name][:], lutData[y_name], OCPI, LUTName)
end

function load_LUT1d!(OCPI::OCPInterface_, filepath::String, LUTName::Symbol; x_name::String="X", y_name::String="Y")
    lutData = matread(filepath)
    OCPI.LUTs1D[LUTName] = LUTs.build1D(lutData[x_name][:], lutData[y_name][:], OCPI, LUTName)
end
