"""
    const OneAtLeast{A, N} = Tuple{A, Vararg{A, N}}

One at least.
"""
const OneAtLeast{A, N} = Tuple{A, Vararg{A, N}}

"""
    const OneOrMore{A} = Union{A, OneAtLeast{A}}

One or more.
"""
const OneOrMore{A} = Union{A, OneAtLeast{A}}

"""
    OneOrMore(x) -> Tuple{typeof(x)}
    OneOrMore(x::Tuple) -> typeof(x)

If `x` is a tuple, return itself; if not, return `(x,)`.
"""
@inline OneOrMore(x) = (x,)
@inline OneOrMore(xs::Tuple) = xs

function ⊞ end
function ⊠ end
function ⊕ end
function ⊗ end
function add! end
function decompose end
function decompose! end
function dimension end
function div! end
function expand end
function expand! end
function id end
function kind end
function permute end
function reset! end
function shape end
function sub! end
function update end
function update! end
function value end
