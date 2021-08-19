module Interfaces

export id, value, rank, dimension
export ⊕, ⊗, ⋅, add!, sub!, mul!, div!
export expand, expand!, decompose, decompose!, permute

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

"Inplace addition."
function add! end

"Inplace subtraction."
function sub! end

"Inplace multiplication."
function mul! end

"Inplace division."
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
