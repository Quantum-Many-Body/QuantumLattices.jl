module AlgebraOverFields

using Printf: @printf
using ..NamedVectors: AbstractNamedVector
using ..TypeTraits: efficientoperations
using ..CompositeStructures: CompositeNTuple,CompositeVector

import ..Utilities.Interfaces: rank,⊗,add!,sub!

export SimpleID,ID
export Element,Elements,idtype
export rank,⊗,add!,sub!

"""
    SimpleID <: AbstractNamedVector

A simple id is the building block of the id system of an algebra over a field.
"""
abstract type SimpleID <: AbstractNamedVector end

"""
    ID(ids::NTuple{N,SimpleID}) where N
    ID(ids::SimpleID...)
    ID(::Type{SID},attrs::Vararg{NTuple{N},M}) where {SID<:SimpleID,N,M}

The id system of an algebra over a field.

Usually, a simple id corresponds to a single generator of the algebra while an id corresponds to an element of the algebra.
"""
struct ID{N,I<:SimpleID} <: CompositeNTuple{N,I}
    contents::NTuple{N,I}
    ID(ids::NTuple{N,SimpleID}) where N=new{N,ids|>eltype}(ids)
end
ID(ids::SimpleID...)=ID(ids)
@generated function ID(::Type{SID},attrs::Vararg{NTuple{N},M}) where {SID<:SimpleID,N,M}
    exprs=[]
    for i=1:N
        args=[:(attrs[$j][$i]) for j=1:M]
        push!(exprs,:(SID($(args...))))
    end
    return :(ID($(exprs...)))
end

"""
    propertynames(::Type{I},private::Bool=false) where I<:ID -> Tuple

Get the property names of a composite id.
"""
Base.propertynames(::Type{I},private::Bool=false) where I<:ID=I|>eltype|>isconcretetype ? cidpropertynames(I,Val(private)) : (:contents,)
@generated function cidpropertynames(::Type{I},::Val{true}) where I<:ID
    @assert I|>eltype|>isconcretetype "cidpropertynames error: not homogeneous input I($(I|>nameof))."
    exprs=[QuoteNode(Symbol(name,'s')) for name in I|>eltype|>fieldnames]
    return Expr(:tuple,QuoteNode(:contents),exprs...)
end
@generated function cidpropertynames(::Type{I},::Val{false}) where I<:ID
    @assert I|>eltype|>isconcretetype "cidpropertynames error: not homogeneous input I($(I|>nameof))."
    exprs=[QuoteNode(Symbol(name,'s')) for name in I|>eltype|>fieldnames]
    return Expr(:tuple,exprs...)
end

"""
    getproperty(cid::ID,name::Symbol)

Get the property of a composite id.
"""
Base.getproperty(cid::ID,name::Symbol)=name==:contents ? getfield(cid,:contents) : cidgetproperty(cid,name)
@generated function cidgetproperty(cid::ID,name::Symbol)
    index=:(index=findfirst(isequal(name),cid|>typeof|>propertynames)::Int)
    exprs=[:(getfield(cid[$i],index)) for i=1:length(cid)]
    return Expr(:block,index,Expr(:tuple,exprs...))
end

"""
    show(io::IO,cid::ID)

Show a composite id.
"""
Base.show(io::IO,cid::ID)=@printf io "%s(%s)" cid|>typeof|>nameof join(cid,",")

"""
    isless(cid1::ID,cid2::ID) -> Bool
    <(cid1::ID,cid2::ID) -> Bool

Compare two ids and judge whether the first is less than the second.

We assume that ids with smaller ranks are always less than those with higher ranks. If two ids are of the same rank, the comparison goes just like that between tuples.
"""
function Base.isless(cid1::ID,cid2::ID)
    r1,r2=cid1|>rank,cid2|>rank
    r1<r2 ? true : r1>r2 ? false : isless(convert(Tuple,cid1),convert(Tuple,cid2))
end
function Base.:<(cid1::ID,cid2::ID)
    r1,r2=cid1|>rank,cid2|>rank
    r1<r2 ? true : r1>r2 ? false : convert(Tuple,cid1)<convert(Tuple,cid2)
end

"""
    rank(::Type{<:ID{N,I}}) where {N,I} -> Int
    rank(id::ID) -> Int

Get the rank of a composite id.
"""
rank(::Type{<:ID{N,I}}) where {N,I}=N
rank(id::ID)=id|>typeof|>rank

"""
    ⊗(sid1::SimpleID,sid2::SimpleID) -> ID
    ⊗(sid::SimpleID,cid::ID) -> ID
    ⊗(cid::ID,sid::SimpleID) -> ID
    ⊗(cid1::ID,cid2::ID) -> ID

Get the direct product of the id system.
"""
⊗(sid1::SimpleID,sid2::SimpleID)=ID(sid1,sid2)
⊗(sid::SimpleID,cid::ID)=ID(sid,convert(Tuple,cid)...)
⊗(cid::ID,sid::SimpleID)=ID(convert(Tuple,cid)...,sid)
⊗(cid1::ID,cid2::ID)=ID(convert(Tuple,cid1)...,convert(Tuple,cid2)...)

"""
    Element{V<:Number,I<:ID}

An element of an algebra over a field.

The first and second attributes of an element must be
- `value::Nuber`: the coefficient of the element
- `id::ID`: the id of the element
"""
abstract type Element{V<:Number,I<:ID,N} end

"""
    valtype(::Type{<:Element{V}}) where {V}
    valtype(m::Element)

Get the type of the value of an element.

The result is also the type of the field over which the algebra is defined.
"""
Base.valtype(::Type{<:Element{V}}) where V=V
Base.valtype(m::Element)=m|>typeof|>valtype

"""
    idtype(::Type{<:Element{V,I}}) where {V,I}
    idtype(m::Element)

The type of the id of an element.
"""
idtype(::Type{<:Element{V,I}}) where {V,I}=I
idtype(m::Element)=m|>typeof|>idtype

"""
    rank(::Type{<:Element{V,I,N}}) where {V,I,N} -> Int
    rank(m::Element) -> Int

Get the rank of an element.
"""
rank(::Type{<:Element{V,I,N}}) where {V,I,N}=N
rank(m::Element)=m|>typeof|>rank

"""
    ==(m1::M,m2::M) where M<:Element -> Bool
    isequal(m1::M,m2::M) where M<:Element -> Bool

Compare two elements and judge whether they are equal to each other.
"""
Base.:(==)(m1::M,m2::M) where M<:Element = ==(efficientoperations,m1,m2)
Base.isequal(m1::M,m2::M) where M<:Element=isequal(efficientoperations,m1,m2)

"""
    replace(m::Element;kwargs...) -> typeof(m)

Return a copy of a concrete `Element` with some of the field values replaced by the keyword arguments.
"""
Base.replace(m::Element;kwargs...)=replace(efficientoperations,m;kwargs...)

"""
    Elements{I<:ID,M<:Element} <: AbstractDict{I,M}

An set of elements of an algebra over a field.

Alias for `Dict{I<:ID,M<:Element}`. Similar iterms are automatically merged thanks to the id system.
"""
const Elements{I<:ID,M<:Element}=Dict{I,M}
"""
    Elements(ms)
    Elements(ms::Pair{I,M}...) where {I<:ID,M<:Element}
    Elements(ms::Element...)

Get the set of elements with similar items merged.
"""
Elements(ms)=Base.dict_with_eltype((K,V)->Dict{K,V},ms,eltype(ms))
function Elements(ms::Pair{I,M}...) where {I<:ID,M<:Element}
    result=Elements{I,M}()
    for (id,m) in ms result[id]=m end
    return result
end
function Elements(ms::Element...)
    result=Elements{ms|>eltype|>idtype,ms|>eltype}()
    for m in ms add!(result,m) end
    return result
end

"""
    zero(ms::Elements) -> typeof(ms)
    zero(::Type{Elements{I,M}}) where {I,M} -> Elements{I,M}

Get a zero set of elements.

A zero set of elements is defined to be the empty one.
"""
Base.zero(ms::Elements)=ms|>typeof|>zero
Base.zero(::Type{Elements{I,M}}) where {I,M}=Elements{I,M}()

"""
    add!(ms::Elements) -> typeof(ms)
    add!(ms::Elements,m::Element) -> typeof(ms)
    add!(ms::Elements,mms::Elements) -> typeof(ms)

Get the inplace addition of elements to a set.
"""
add!(ms::Elements)=ms
function add!(ms::Elements,m::Element)
    @assert ms|>valtype==m|>typeof "add! error: dismatched type, $(ms|>valtype) and $(m|>typeof)."
    mid=m.id
    ms[mid]=haskey(ms,mid) ? replace(m,value=ms[mid].value+m.value) : m
    abs(ms[mid].value)==0.0 && delete!(ms,mid)
    ms
end
add!(ms::Elements,mms::Elements)=(for m in mms|>values add!(ms,m) end; ms)

"""
    sub!(ms::Elements) -> typeof(ms) -> typeof(ms)
    sub!(ms::Elements,m::Element) -> typeof(ms)
    sub!(ms::Elements,mms::Elements) -> typeof(ms)

Get the inplace subtraction of elements from a set.
"""
sub!(ms::Elements)=ms
function sub!(ms::Elements,m::Element)
    @assert ms|>valtype==m|>typeof "sub! error: dismatched type, $(ms|>valtype) and $(m|>typeof)."
    mid=m.id
    ms[mid]=haskey(ms,mid) ? replace(m,value=ms[mid].value-m.value) : -m
    abs(ms[mid].value)==0.0 && delete!(ms,mid)
    ms
end
sub!(ms::Elements,mms::Elements)=(for m in mms|>values sub!(ms,m) end; ms)

"""
    +(m::Element) -> typeof(m)
    +(ms::Elements) -> typeof(ms)
    +(ms::Elements,m::Element) -> Elements
    +(m1::Element,m2::Element) -> Elements
    +(m::Element,ms::Elements) -> Elements
    +(ms1::Elements,ms2::Elements) -> Elements

Overloaded `+` operator between elements of an algebra over a field.
"""
Base.:+(m::Element)=m
Base.:+(ms::Elements)=ms
Base.:+(ms::Elements,m::Element)=m+ms
Base.:+(m1::Element,m2::Element)=add!(Elements{typejoin(m1|>idtype,m2|>idtype),typejoin(m1|>typeof,m2|>typeof)}(m1.id=>m1),m2)
Base.:+(m::Element,ms::Elements)=add!(Elements{typejoin(m|>idtype,ms|>keytype),typejoin(m|>typeof,ms|>valtype)}(ms),m)
Base.:+(ms1::Elements,ms2::Elements)=add!(Elements{typejoin(ms1|>keytype,ms2|>keytype),typejoin(ms1|>valtype,ms2|>valtype)}(ms1),ms2)

"""
    *(factor::Number,m::Element) -> Element
    *(m::Element,factor::Number) -> Element
    *(m1::Element,m2::Element) -> Element
    *(factor::Number,ms::Elements) -> Elements
    *(ms::Elements,factor::Number) -> Elements
    *(m::Element,ms::Elements) -> Elements
    *(ms::Elements,m::Element) -> Elements
    *(ms1::Elements,ms2::Elements) -> Elements

Overloaded `*` operator for element-scalar multiplications and element-element multiplications of an algebra over a field.
"""
Base.:*(factor::Number,m::Element)=m*factor
Base.:*(factor::Number,ms::Elements)=ms*factor
Base.:*(m::Element,factor::Number)=replace(m,value=m.value*factor)
Base.:*(ms::Elements,factor::Number)=abs(factor)==0.0 ? zero(Elements) : Elements(id=>m*factor for (id,m) in ms)
Base.:*(m::Element,ms::Elements)=Elements((m*mm for mm in ms|>values)...)
Base.:*(ms::Elements,m::Element)=Elements((mm*m for mm in ms|>values)...)
Base.:*(ms1::Elements,ms2::Elements)=Elements((m1*m2 for m1 in ms1|>values for m2 in ms2|>values)...)
function Base.:*(m1::Element,m2::Element)
    @assert(    m1|>typeof|>nameof==m2|>typeof|>nameof && m1|>typeof|>fieldcount==m2|>typeof|>fieldcount==2,
                "\"*\" error: not implemented between $(m1|>typeof|>nameof) and $(m2|>typeof|>nameof)."
                )
    typeof(m1).name.wrapper(m1.value*m2.value,m1.id⊗m2.id)
end

"""
    -(m::Element) -> typeof(m)
    -(ms::Elements) -> typeof(ms)
    -(m1::Element,m2::Element) -> Elements
    -(m::Element,ms::Elements) -> Elements
    -(ms::Elements,m::Element) -> Elements
    -(ms1::Elements,ms2::Elements) -> Elements

Overloaded `-` operator between elements of an algebra over a field.
"""
Base.:-(m::Element)=m*(-1)
Base.:-(ms::Elements)=ms*(-1)
Base.:-(m1::Element,m2::Element)=sub!(Elements{typejoin(m1|>idtype,m2|>idtype),typejoin(m1|>typeof,m2|>typeof)}(m1.id=>m1),m2)
Base.:-(m::Element,ms::Elements)=sub!(Elements{typejoin(m|>idtype,ms|>keytype),typejoin(m|>typeof,ms|>valtype)}(m.id=>m),ms)
Base.:-(ms::Elements,m::Element)=sub!(Elements{typejoin(m|>idtype,ms|>keytype),typejoin(m|>typeof,ms|>valtype)}(ms),m)
Base.:-(ms1::Elements,ms2::Elements)=sub!(Elements{typejoin(ms1|>keytype,ms2|>keytype),typejoin(ms1|>valtype,ms2|>valtype)}(ms1),ms2)

"""
    /(m::Element,factor::Number)
    /(ms::Elements,factor::Number)

Overloaded `/` operator for element-sclar division of an algebra over a field.
"""
Base.:/(m::Element,factor::Number)=m*(1/factor)
Base.:/(ms::Elements,factor::Number)=ms*(1/factor)

end #module
