module build

using ModelingToolkit, JuMP
using LinearAlgebra
using DataInterpolations
using DataInterpolationsND
import JuMP: add_nonlinear_operator

include("LGPoly.jl")
include("LUTs.jl")

include("build/LUTs/types.jl")
include("build/types.jl")

include("build/architecture.jl")
include("build/settings.jl")
include("build/bounds.jl")
include("build/variables.jl")
include("build/expression.jl")
include("build/constraints.jl")
include("build/results.jl")


include("build/LUTs/builders.jl")
include("build/LUTs/loader.jl")

end
