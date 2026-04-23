function Interpolation2d(; itp, name)
    @parameters (interpolator2d::typeof(itp))(..) = itp
    @named input1 = RealInput()
    @named input2 = RealInput()
    @named output = RealOutput()

    eqs = [output.u ~ interpolator2d(input1.u, input2.u)]

    return System(eqs, t, [], [interpolator2d]; name, systems = [input1,input2, output])
end

function Interpolation1d(; itp, name)
    @parameters (interpolator1d::typeof(itp))(..) = itp
    @named input = RealInput()
    @named output = RealOutput()

    eqs = [output.u ~ interpolator1d(input.u)]

    return System(eqs, t, [], [interpolator1d]; name, systems = [input, output])
end
