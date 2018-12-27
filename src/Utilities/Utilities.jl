module Utilities

using Formatting: FormatSpec,fmt

export atol,rtol,Float
export forder,corder,indtosub,subtoind
export decimaltostr,ordinal,efficientoperations,delta

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

@generated tail(ts::NTuple{N}) where N=Expr(:tuple,[:(ts[$i]) for i=2:N]...)
@generated head(ts::NTuple{N}) where N=Expr(:tuple,[:(ts[$i]) for i=1:N-1]...)

"""
    indtosub(dims::Tuple,ind::Int,order::FOrder) -> Tuple
    indtosub(dims::Tuple,ind::Int,order::COrder) -> Tuple

Convert an linear index to Cartesian index. Fortran-order or C-order can be assigned.
"""
function indtosub(dims::Tuple,ind::Int,order::FOrder)
    length(dims)==0 && return ()
    ((ind-1)%dims[1]+1,indtosub(tail(dims),(ind-1)÷dims[1]+1,order)...)
end
function indtosub(dims::Tuple,ind::Int,order::COrder)
    length(dims)==0 && return ()
    (indtosub(head(dims),(ind-1)÷dims[end]+1,order)...,(ind-1)%dims[end]+1)
end

"""
    subtoind(dims::NTuple{N,Int},inds::NTuple{N,Int},order::FOrder) where N -> Int
    subtoind(dims::NTuple{N,Int},inds::NTuple{N,Int},order::COrder) where N -> Int

Convert an Cartesian index to linear index. Fortran-order or C-order can be assigned.
"""
function subtoind(dims::NTuple{N,Int},inds::NTuple{N,Int},order::FOrder) where N
    length(dims)==0 && return 1
    (subtoind(tail(dims),tail(inds),order)-1)*dims[1]+inds[1]
end
function subtoind(dims::NTuple{N,Int},inds::NTuple{N,Int},order::COrder) where N
    length(dims)==0 && return 1
    (subtoind(head(dims),head(inds),order)-1)*dims[end]+inds[end]
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

struct EfficientOperations end
"""
    efficientoperations

Indicate that the efficient operations, i.e. "=="/"isequal", "<"/"isless" or "replace", will be used.
"""
const efficientoperations=EfficientOperations()

"""
    ==(::EfficientOperations,o1,o2) -> Bool
    isequal(::EfficientOperations,o1,o2) -> Bool

Compare two objects and judge whether they are eqaul to each other.
"""
@generated function Base.:(==)(::EfficientOperations,o1,o2)
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
@generated function Base.isequal(::EfficientOperations,o1,o2)
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
    <(::EfficientOperations,o1,o2) -> Bool
    isless(::EfficientOperations,o1,o2) -> Bool

Compare two objects and judge whether the first is less than the second.
"""
@generated function Base.:<(::EfficientOperations,o1,o2)
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
@generated function Base.isless(::EfficientOperations,o1,o2)
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

"""
    replace(::EfficientOperations,o;kwargs...) -> typeof(o)

Return a copy of the input object with some of the field values replaced by the keyword arguments.
"""
@generated function Base.replace(::EfficientOperations,o;kwargs...)
    exprs=[:(get(kwargs,$name,getfield(o,$name))) for name in QuoteNode.(o|>fieldnames)]
    return :(typeof(o).name.wrapper($(exprs...)))
end

"""
    delta(i,j) -> Int

Kronecker delta function.
"""
delta(i,j)=i==j ? 1 : 0

"Generic interface of the direct sum of some types."
function ⊕ end

"Generic interface of the direct product of some types."
function ⊗ end

"Generic interface of the inplace addition of some types."
function add! end

"Generic interface of the inplace subtraction of some types."
function sub! end

"Generic interface of the inplace multiplication of some types."
function mul! end

"Generic interface of the inplace division of some types."
function div! end

"Generic interface of the rank of some types."
function rank end

"Generic interface of the dimension of some types."
function dimension end

"Generic interface of the expansion of some types."
function expand end

"Generic interface of the permutation of some types."
function permute end

"Generic interface of the vector representation of some types."
function vector end

"Generic interface of the matrix representation of some types."
function matrix end

include("Factory.jl")
include("CompositeStructure.jl")
include("Tree.jl")
include("NamedVector.jl")
include("AlgebraOverField.jl")
include("QuantumNumber.jl")

end # module
