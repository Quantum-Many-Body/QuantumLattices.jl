module Essentials

using Reexport: @reexport

include("Spatial.jl")
include("DegreeOfFreedom.jl")

@reexport using .Spatial
@reexport using .DegreeOfFreedom

end # module
