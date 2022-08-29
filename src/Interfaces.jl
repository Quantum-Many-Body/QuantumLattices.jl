module Interfaces

import LinearAlgebra: rank, mul!, ⋅

export dimension, id, rank, value
export ⊕, ⊗, ⋅, add!, div!, mul!, sub!
export decompose, decompose!, expand, expand!, permute

"""
ID.
"""
function id end

"""
Value.
"""
function value end

"Rank."
function rank end

"Dimension."
function dimension end

"Direct sum."
function ⊕ end

"Direct product."
function ⊗ end

"Dot product."
function ⋅ end

"In-place addition."
function add! end

"In-place subtraction."
function sub! end

"In-place multiplication."
function mul! end

"In-place division."
function div! end

"Get the expansion."
function expand end

"In place expansion."
function expand! end

"Decompose."
function decompose end

"In place decomposition."
function decompose! end

"Get the permutation."
function permute end

end # module
