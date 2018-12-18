module Utilities

using Formatting: FormatSpec,fmt

export atol,rtol,Float
export forder,corder,ind2sub,sub2ind
export decimaltostr,ordinal,comparison

"Absolute tolerance for float numbers."
const atol=5*10^-14

"Relative tolerance for float numbers."
const rtol=√atol

"Default float type."
const Float=Float64

abstract type MemoryOrder end

struct FOrder <: MemoryOrder end
"""
    forder

Indicate that the convertion between Cartesian index and linear index is using the Fortran order.
"""
const forder=FOrder()

struct COrder <: MemoryOrder end
"""
    corder

Indicate that the convertion between Cartesian index and linear index is using the C/C++ order.
"""
const corder=COrder()

"""
    ind2sub(dims::Tuple,ind::Int,order::FOrder) -> Tuple
    ind2sub(dims::Tuple,ind::Int,order::COrder) -> Tuple

Convert an linear index to Cartesian index. Fortran-order or C-order can be assigned.
"""
function ind2sub(dims::Tuple,ind::Int,order::FOrder)
    length(dims)==0 && return ()
    ((ind-1)%dims[1]+1,ind2sub(dims[2:end],(ind-1)÷dims[1]+1,order)...)
end
function ind2sub(dims::Tuple,ind::Int,order::COrder)
    length(dims)==0 && return ()
    (ind2sub(dims[1:end-1],(ind-1)÷dims[end]+1,order)...,(ind-1)%dims[end]+1)
end

"""
    sub2ind(dims::NTuple{N,Int},inds::NTuple{N,Int},order::FOrder) where N -> Int
    sub2ind(dims::NTuple{N,Int},inds::NTuple{N,Int},order::COrder) where N -> Int

Convert an Cartesian index to linear index. Fortran-order or C-order can be assigned.
"""
function sub2ind(dims::NTuple{N,Int},inds::NTuple{N,Int},order::FOrder) where N
    length(dims)==0 && return 1
    (sub2ind(dims[2:end],inds[2:end],order)-1)*dims[1]+inds[1]
end
function sub2ind(dims::NTuple{N,Int},inds::NTuple{N,Int},order::COrder) where N
    length(dims)==0 && return 1
    (sub2ind(dims[1:end-1],inds[1:end-1],order)-1)*dims[end]+inds[end]
end

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

struct Comparison end
"""
    comparison

Indicate that the comparison methods, i.e. "=="/"isequal" or "<"/"isless", will be used.
"""
const comparison=Comparison()

"""
    ==(::Comparison,o1,o2) -> Bool
    isequal(::Comparison,o1,o2) -> Bool

Compare two objects and judge whether they are eqaul to each other.
"""
@generated function Base.:(==)(::Comparison,o1,o2)
    fcount=o1|>fieldcount
    if fcount==o2|>fieldcount
        if fcount==0
            return :(true)
        else
            expr=:(getfield(o1,1)==getfield(o2,1))
            for i=2:fcount
                expr=Expr(:&&,expr,:(getfield(o1,$i)==getfield(o2,$i)))
            end
            return expr
        end
    else
        return :(false)
    end
end
@generated function Base.isequal(::Comparison,o1,o2)
    fcount=o1|>fieldcount
    if fcount==o2|>fieldcount
        if fcount==0
            return :(true)
        else
            expr=:(isequal(getfield(o1,1),getfield(o2,1)))
            for i=2:fcount
                expr=Expr(:&&,expr,:(isequal(getfield(o1,$i),getfield(o2,$i))))
            end
            return expr
        end
    else
        return :(false)
    end
end

"""
    <(::Comparison,o1,o2) -> Bool
    isless(::Comparison,o1,o2) -> Bool

Compare two objects and judge whether the first is less than the second.
"""
@generated function Base.:<(::Comparison,o1,o2)
    n1,n2=o1|>fieldcount,o2|>fieldcount
    N=min(n1,n2)
    expr=n1<n2 ? Expr(:if,:(getfield(o1,$N)==getfield(o2,$N)),true,false) : false
    expr=Expr(:if,:(getfield(o1,$N)<getfield(o2,$N)),true,expr)
    for i in range(N-1,stop=1,step=-1)
        expr=Expr(:if,:(getfield(o1,$i)>getfield(o2,$i)),false,expr)
        expr=Expr(:if,:(getfield(o1,$i)<getfield(o2,$i)),true,expr)
    end
    return expr
end
@generated function Base.isless(::Comparison,o1,o2)
    n1,n2=o1|>fieldcount,o2|>fieldcount
    N=min(n1,n2)
    expr=n1<n2 ? Expr(:if,:(getfield(o1,$N)==getfield(o2,$N)),true,false) : false
    expr=Expr(:if,:(isless(getfield(o1,$N),getfield(o2,$N))),true,expr)
    for i in range(N-1,stop=1,step=-1)
        expr=Expr(:if,:(isless(getfield(o2,$i),getfield(o1,$i))),false,expr)
        expr=Expr(:if,:(isless(getfield(o1,$i),getfield(o2,$i))),true,expr)
    end
    return expr
end

"Generic interface of the direct sum of some types."
function ⊕ end

"Generic interface of the direct product of some types."
function ⊗ end

"Generic interface of the rank of some types."
function rank end

"Generic interface of the dimension of some types."
function dimension end

"Generic interface of permuting of some types."
function permute end

include("Factory.jl")
include("CompositeStructure.jl")
include("Tree.jl")
include("NamedVector.jl")
include("AlgebraOverField.jl")
include("QuantumNumber.jl")

end # module
