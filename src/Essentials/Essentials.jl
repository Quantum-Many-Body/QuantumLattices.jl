module Essentials

using Reexport: @reexport

export dtype, kind, update!, reset!

function dtype end
function kind end
function update! end
function reset! end

include("QuantumAlgebras.jl")
include("QuantumNumbers.jl")
include("Spatials.jl")
include("DegreesOfFreedom.jl")
include("Terms.jl")
include("Frameworks.jl")
include("QuantumSystems.jl")

@reexport using .QuantumAlgebras
@reexport using .QuantumNumbers
@reexport using .Spatials
@reexport using .DegreesOfFreedom
@reexport using .Terms
@reexport using .Frameworks
@reexport using .QuantumSystems

end # module
