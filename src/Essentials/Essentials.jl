module Essentials

using Reexport: @reexport

export dtype, kind, update!, reset!

function dtype end
function kind end
function update! end
function reset! end

include("Spatials.jl")
include("DegreesOfFreedom.jl")
include("Terms.jl")
include("Frameworks.jl")
include("FockPackage.jl")
include("SpinPackage.jl")
include("PhononPackage.jl")

@reexport using .Spatials
@reexport using .DegreesOfFreedom
@reexport using .Terms
@reexport using .Frameworks
@reexport using .FockPackage
@reexport using .SpinPackage
@reexport using .PhononPackage

end # module
