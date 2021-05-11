module Essentials

using Reexport: @reexport

export kind, update!, reset!

function kind end
function update! end
function reset! end

include("Spatials.jl")
include("DegreesOfFreedom.jl")
include("Terms.jl")
include("Frameworks.jl")
include("FockPackage.jl")
include("SpinPackage.jl")

@reexport using .Spatials
@reexport using .DegreesOfFreedom
@reexport using .Terms
@reexport using .Frameworks
@reexport using .FockPackage
@reexport using .SpinPackage

end # module
