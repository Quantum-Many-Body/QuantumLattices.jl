module Essentials

using Reexport: @reexport

include("Spatial.jl")

@reexport using .Spatial

end # module
