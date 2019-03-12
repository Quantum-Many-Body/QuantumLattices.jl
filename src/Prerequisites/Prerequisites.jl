module Prerequisites

using Formatting: FormatSpec,fmt

export atol,rtol,Float
export decimaltostr,ordinal,delta

"Absolute tolerance for float numbers."
const atol=5*10^-14

"Relative tolerance for float numbers."
const rtol=√atol

"Default float type."
const Float=Float64

"""
    decimaltostr(number::Integer,n::Int=5)
    decimaltostr(number::Rational,n::Int=5)
    decimaltostr(number::AbstractFloat,n::Int=5)
    decimaltostr(number::Complex,n::Int=5)

Convert a number to a string with at most `n` decimal places.
"""
decimaltostr(number::Integer,::Int=5)=string(number)
decimaltostr(number::Rational,::Int=5)=string(number)
function decimaltostr(number::AbstractFloat,n::Int=5)
    if number==0.0
        result="0.0"
    elseif 10^-5<abs(number)<10^6
        result=rstrip(fmt(FormatSpec(".$(n)f"),number),'0')
        result[end]=='.' && (result=result*'0')
    else
        result=fmt(FormatSpec(".$(n)e"),number)
        epos=findfirst(isequal('e'),result)
        temp=rstrip(result[1:epos-1],'0')
        result=temp[end]=='.' ? temp*"0"*result[epos:end] : temp*result[epos:end]
    end
    result
end
function decimaltostr(number::Complex,n::Int=5)
    sreal=real(number)==0 ? "0" : decimaltostr(real(number),n)
    simag=imag(number)==0 ? "0" : decimaltostr(imag(number),n)
    result=""
    sreal=="0" || (result=result*sreal)
    simag=="0" || (result=(simag[1]=='-' ? result*simag : length(result)==0 ? simag : result*"+"*simag)*"im")
    length(result)==0 && (result="0.0")
    return result
end

"""
    ordinal(number::Interger)

Convert a positive number to its corresponding ordinal.
"""
function ordinal(number)
    @assert number>0
    number==1 ? "1st" : number==2 ? "2nd" : number==3 ? "3rd" : "$(number)th"
end

"""
    delta(i,j) -> Int

Kronecker delta function.
"""
delta(i,j)=i==j ? 1 : 0

include("TypeTraits.jl")
include("Factories.jl")
include("CompositeStructures.jl")
include("Trees.jl")
include("NamedVectors.jl")

end # module
