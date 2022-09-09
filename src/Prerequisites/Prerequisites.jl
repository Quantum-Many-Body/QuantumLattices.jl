module Prerequisites

using Formatting: FormatSpec, fmt
using Reexport: @reexport

export atol, rtol, Float
export allequal, concatenate, decimaltostr, delta, ordinal

"Absolute tolerance for float numbers."
const atol = 5 * 10^-14

"Relative tolerance for float numbers."
const rtol = √atol

"Default float type."
const Float = Float64

"""
    allequal(xs) -> Bool

Judge whether all elements are equal.
"""
@inline allequal(xs) = all(isequal(first(xs)), xs)

"""
    decimaltostr(number, ::Int=5)
    decimaltostr(number::Integer, n::Int=5)
    decimaltostr(number::Rational, n::Int=5)
    decimaltostr(number::AbstractFloat, n::Int=5)
    decimaltostr(number::Complex, n::Int=5)

Convert a number to a string with at most `n` decimal places.
"""
@inline decimaltostr(number, ::Int=5) = repr(number)
@inline decimaltostr(number::Integer, ::Int=5) = string(number)
@inline decimaltostr(number::Rational, ::Int=5) = string(number)
function decimaltostr(number::AbstractFloat, n::Int=5)
    if number == 0.0
        result = "0.0"
    elseif 10^-5 < abs(number) < 10^6
        result = rstrip(fmt(FormatSpec(".$(n)f"), number), '0')
        (result[end] == '.') && (result = result * '0')
    else
        result = fmt(FormatSpec(".$(n)e"), number)
        epos = findfirst(isequal('e'), result)
        temp = rstrip(result[1:epos-1], '0')
        result = (temp[end] == '.') ? (temp * "0" * result[epos:end]) : (temp * result[epos:end])
    end
    result
end
function decimaltostr(number::Complex, n::Int=5)
    sreal = (real(number) == 0) ? "0" : decimaltostr(real(number), n)
    simag = (imag(number) == 0) ? "0" : decimaltostr(imag(number), n)
    result = ""
    (sreal == "0") || (result = result * sreal)
    (simag == "0") || (result = ((simag[1] == '-') ? (result * simag) : (length(result) == 0) ? simag : (result * "+" * simag)) * "im")
    (length(result) == 0) && (result = "0.0")
    return result
end

"""
    ordinal(number::Integer)

Convert a positive number to its corresponding ordinal.
"""
@inline function ordinal(number::Integer)
    @assert number > 0 "ordinal error: input number must be positive."
    (number == 1) ? "1st" : (number == 2) ? "2nd" : (number == 3) ? "3rd" : "$(number)th"
end

"""
    delta(i, j) -> Int

Kronecker delta function.
"""
@inline delta(i, j) = (i == j) ? 1 : 0

"""
    concatenate(ts::Tuple...) -> Tuple

Concatenate tuples.
"""
@inline @generated concatenate(ts::Tuple...) = Expr(:tuple, map(i->:(ts[$i]...), 1:length(ts))...)

"""
    searchsortedfirst(table, basis; compare=<) -> Int

Use the binary search method to find the position of a basis in a sorted table so that the order is preserved if the basis in inserted in that position.
"""
function Base.searchsortedfirst(table, basis; compare=<)
    lo, hi = 0, length(table)+1
    @inbounds while lo < hi-1
        m = (lo+hi) >>> 1
        if compare(table[m], basis)
            lo = m
        else
            hi = m
        end
    end
    return hi
end

include("Combinatorics.jl")
include("Traits.jl")
include("CompositeStructures.jl")
include("SimpleTrees.jl")
include("NamedVectors.jl")
include("VectorSpaces.jl")

@reexport using .Combinatorics
@reexport using .Traits
@reexport using .CompositeStructures
@reexport using .SimpleTrees
@reexport using .NamedVectors
@reexport using .VectorSpaces

end # module
