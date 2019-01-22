module Essentials

using Reexport: @reexport

include("Spatials.jl")
include("DegreesOfFreedom.jl")
include("FockPackage/FockPackage.jl")
include("SpinPackage/SpinPackage.jl")

@reexport using .Spatials
@reexport using .DegreesOfFreedom
@reexport using .FockPackage
@reexport using .SpinPackage

end # module
