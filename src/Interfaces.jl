module Interfaces

export ⊕,⊗,⋅
export add!,sub!,mul!,div!
export id,rank,dimension,kind,sequence,dims,index,inds
export update!,reset!,register!,prepare!,run!,expand,expand!,decompose,decompose!,regularize,regularize!
export permute,vector,matrix

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

"Inplace Update."
function update! end

"Inplace reset."
function reset! end

"Inplace register."
function register! end

"Inplace prepare."
function prepare! end

"Inplace run."
function run! end

"Id."
function id end

"Rank."
function rank end

"Dimension."
function dimension end

"Kind."
function kind end

"Sequence."
function sequence end

"Dimensions."
function dims end

"Index."
function index end

"Indices."
function inds end

"Get the expansion."
function expand end

"In place expansion."
function expand! end

"Decompose."
function decompose end

"In place decomposition."
function decompose! end

"Regularize."
function regularize end

"In place regularization."
function regularize! end

"Get the permutation."
function permute end

"Vector representation."
function vector end

"Matrix representation."
function matrix end

end # module
