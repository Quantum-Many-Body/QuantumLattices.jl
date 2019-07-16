module AlgebraOverFields

using Printf: @printf,@sprintf
using ...Interfaces: dimension
using ...Prerequisites: atol,rtol,Float,rawtype
using ...Prerequisites.NamedVectors: NamedVector
using ...Prerequisites.TypeTraits: efficientoperations
using ..Combinatorics: AbstractCombinatorics

import ...Interfaces: rank,add!,sub!,mul!,div!,⊗,⋅,sequence,permute

export SimpleID,ID
export idpropertynames,idisless
export IdSpace
export Element,ScalarElement,Elements
export scalartype,idtype

"""
    SimpleID <: NamedVector

A simple id is the building block of the id system of an algebra over a field.
"""
abstract type SimpleID <: NamedVector end

"""
    ID{I<:SimpleID,N}

The (composite) id system of an algebra over a field.

Type alias for `NTuple{N,I} where {N,I<:SimpleID}`
"""
const ID{I<:SimpleID,N}=NTuple{N,I}
"""
    ID(ids::SimpleID...)
    ID(ids::NTuple{N,SimpleID}) where N

Get the composite id from simple ids.
"""
ID(ids::SimpleID...)=ids
ID(ids::NTuple{N,SimpleID}) where N=ids
"""
    ID(::Type{SID},attrs::Vararg{NTuple{N},M}) where {SID<:SimpleID,N,M}

Get the composite id from the components of simple ids.
"""
@generated function ID(::Type{SID},attrs::Vararg{NTuple{N,Any},M}) where {SID<:SimpleID,N,M}
    exprs=[]
    for i=1:N
        args=[:(attrs[$j][$i]) for j=1:M]
        push!(exprs,:(SID($(args...))))
    end
    return :(ID($(exprs...)))
end

"""
    idpropertynames(::Type{I}) where I<:ID{SimpleID} -> Tuple{Vararg{Symbol}}

Get the property names of a composite id.
"""
@generated function idpropertynames(::Type{I}) where I<:ID{SimpleID}
    exprs=[QuoteNode(Symbol(name,'s')) for name in I|>eltype|>fieldnames]
    return Expr(:tuple,exprs...)
end

"""
    getproperty(cid::ID{SimpleID},name::Symbol)

Get the property of a composite id.
"""
Base.getproperty(cid::ID{SimpleID},name::Symbol)=idgetproperty(cid,Val(name),Val(cid|>typeof|>idpropertynames))
@generated function idgetproperty(cid::ID{SimpleID},::Val{name},::Val{names}) where {name,names}
    index=findfirst(isequal(name),names)::Int
    exprs=[:(getfield(cid[$i],$index)) for i=1:fieldcount(cid)]
    return Expr(:tuple,exprs...)
end

"""
    show(io::IO,cid::Tuple{SimpleID,Vararg{SimpleID}})

Show a composite id.
"""
Base.show(io::IO,cid::Tuple{SimpleID,Vararg{SimpleID}})=@printf io "ID(%s)" join(cid,",")

"""
    promote_type(::Type{Tuple{}},::Type{Tuple{}})
    promote_type(::Type{Tuple{}},I::Type{<:ID{SimpleID,N}}) where N
    promote_type(I::Type{<:ID{SimpleID,N}},::Type{Tuple{}}) where N

Define the promote rule for ID types.
"""
Base.promote_type(::Type{Tuple{}},::Type{Tuple{}})=Tuple{}
Base.promote_type(::Type{Tuple{}},I::Type{<:ID{SimpleID,N}}) where N=ID{I|>eltype}
Base.promote_type(I::Type{<:ID{SimpleID,N}},::Type{Tuple{}}) where N=ID{I|>eltype}

"""
    idisless(cid1::ID{SimpleID},cid2::ID{SimpleID}) -> Bool

Compare two ids and judge whether the first is less than the second.

We define that ids with smaller ranks are always less than those with higher ranks. If two ids are of the same rank, the comparison goes just like that between tuples.
"""
function idisless(cid1::ID{SimpleID},cid2::ID{SimpleID})
    r1,r2=cid1|>rank,cid2|>rank
    r1<r2 ? true : r1>r2 ? false : isless(cid1,cid2)
end

"""
    rank(id::ID{SimpleID}) -> Int
    rank(::Type{<:ID{SimpleID}}) -> Any
    rank(::Type{<:ID{SimpleID,N}}) where N -> Int

Get the rank of a composite id.
"""
rank(id::ID{SimpleID})=id|>typeof|>rank
rank(::Type{<:ID{SimpleID}})=Any
rank(::Type{<:ID{SimpleID,N}}) where N=N

"""
    *(sid1::SimpleID,sid2::SimpleID) -> ID{SimpleID}
    *(sid::SimpleID,cid::ID{SimpleID}) -> ID{SimpleID}
    *(cid::ID{SimpleID},sid::SimpleID) -> ID{SimpleID}
    *(cid1::ID{SimpleID},cid2::ID{SimpleID}) -> ID{SimpleID}

Get the product of the id system.
"""
Base.:*(sid1::SimpleID,sid2::SimpleID)=ID(sid1,sid2)
Base.:*(sid::SimpleID,cid::ID{SimpleID})=ID(sid,cid...)
Base.:*(cid::ID{SimpleID},sid::SimpleID)=ID(cid...,sid)
Base.:*(cid1::ID{SimpleID},cid2::ID{SimpleID})=ID(cid1...,cid2...)

"""
    Element{V,I<:ID{SimpleID}}

An element of an algebra over a field.

The first and second attributes of an element must be
- `value`: the coefficient of the element
- `id::ID{SimpleID}`: the id of the element
"""
abstract type Element{V,I<:ID{SimpleID}} end

"""
    ScalarElement{V}

Identity element.
"""
const ScalarElement{V}=Element{V,Tuple{}}
scalartype(::Type{M}) where M<:Element=rawtype(M){M|>valtype,Tuple{}}

"""
    valtype(m::Element)
    valtype(::Type{<:Element})
    valtype(::Type{<:Element{V}}) where V

Get the type of the value of an element.

The result is also the type of the field over which the algebra is defined.
"""
Base.valtype(m::Element)=m|>typeof|>valtype
Base.valtype(::Type{<:Element})=Any
Base.valtype(::Type{<:Element{V}}) where V=V

"""
    idtype(m::Element)
    idtype(::Type{<:Element{V,I} where V}) where I<:ID{SimpleID}
    idtype(::Type{<:Element})

The type of the id of an element.
"""
idtype(m::Element)=m|>typeof|>idtype
idtype(::Type{<:Element{V,I} where V}) where I<:ID{SimpleID}=I
idtype(::Type{<:Element})=ID{SimpleID}

"""
    rank(m::Element) -> Int
    rank(::Type{M}) where M<:Element -> Int

Get the rank of an element.
"""
rank(m::Element)=m|>typeof|>rank
rank(::Type{M}) where M<:Element=M|>idtype|>rank

"""
    ==(m1::Element,m2::Element) -> Bool
    isequal(m1::Element,m2::Element) -> Bool

Compare two elements and judge whether they are equal to each other.
"""
Base.:(==)(m1::Element,m2::Element) = ==(efficientoperations,m1,m2)
Base.isequal(m1::Element,m2::Element)=isequal(efficientoperations,m1,m2)

"""
    isapprox(m1::Element,m2::Element;atol::Real=atol,rtol::Real=rtol) -> Bool

Compare two elements and judge whether they are inexactly equivalent to each other.
"""
Base.isapprox(m1::Element,m2::Element;atol::Real=atol,rtol::Real=rtol)=isapprox(efficientoperations,Val((:value,)),m1,m2;atol=atol,rtol=rtol)

"""
    replace(m::Element;kwargs...) -> typeof(m)

Return a copy of a concrete `Element` with some of the field values replaced by the keyword arguments.
"""
Base.replace(m::Element;kwargs...)=replace(efficientoperations,m;kwargs...)

"""
    promote_rule(::Type{M1},::Type{M2}) where {M1<:Element,M2<:Element}

Define the promote rule for Element types.
"""
function Base.promote_rule(::Type{M1},::Type{M2}) where {M1<:Element,M2<:Element}
    M1<:M2 && return M2
    M2<:M1 && return M1
    v1,i1,r1=M1|>valtype,M1|>idtype,M1|>rank
    v2,i2,r2=M2|>valtype,M2|>idtype,M2|>rank
    V,I=promote_type(v1,v2),promote_type(i1,i2)
    M=r2==0 ? rawtype(M1) : r1==0 ? rawtype(M2) : typejoin(rawtype(M1),rawtype(M2))
    return isconcretetype(I) ? M{V,I} : M{V,<:I}
end

"""
    mul_promote_type(::Type{M},::Type{V}) where {M<:Element,V<:Number}
    mul_promote_type(::Type{V},::Type{M}) where {M<:Element,V<:Number}

Define the promote rule for the multiplication between an Element and a scalar.
"""
mul_promote_type(::Type{M},::Type{V}) where {M<:Element,V<:Number}=mul_promote_type(V,M)
function mul_promote_type(::Type{V},::Type{M}) where {M<:Element,V<:Number}
    C,I=promote_type(M|>valtype,V),M|>idtype
    return isconcretetype(I) ? rawtype(M){C,I} : rawtype(M){C,<:I}
end

"""
    one(::Type{M}) where {M<:Element}

Get the identity operator.
"""
function Base.one(::Type{M}) where {M<:Element}
    rtype=rawtype(M)
    vtype=isconcretetype(valtype(M)) ? valtype(M) : Int
    @assert fieldnames(rtype)==(:value,:id) "one error: not supproted type($(nameof(rtype)))."
    return rtype(one(vtype),ID())
end

"""
    convert(::Type{M},m::ScalarElement) where {M<:ScalarElement}
    convert(::Type{M},m::Number) where {M<:ScalarElement}
    convert(::Type{M},m::Element) where {M<:Element}

1) Convert a scalar element from one type to another;
2) Convert a scalar to a scalar element;
3) Convert an element from one type to another.
"""
function Base.convert(::Type{M},m::ScalarElement) where {M<:ScalarElement}
    typeof(m)<:M && return m
    @assert fieldnames(M)==(:value,:id) "convert error: not supported type($(nameof(M)))."
    return rawtype(M)(convert(M|>valtype,m.value),m.id)
end
function Base.convert(::Type{M},m::Number) where {M<:ScalarElement}
    @assert fieldnames(M)==(:value,:id) "convert error: not supported type($(nameof(M)))."
    return rawtype(M)(convert(M|>valtype,m),ID())
end
function Base.convert(::Type{M},m::Element) where {M<:Element}
    typeof(m)<:M && return m
    @assert idtype(m)<:idtype(M) "convert error: dismatched ID type."
    @assert rawtype(typeof(m))<:rawtype(M) "convert error: dismatched raw Element type."
    return replace(m,value=convert(valtype(M),m.value))
end

"""
    Elements{I<:ID{SimpleID},M<:Element} <: AbstractDict{I,M}

An set of elements of an algebra over a field.

Type alias for `Dict{I<:ID{SimpleID},M<:Element}`.
Similar iterms are automatically merged thanks to the id system.
"""
const Elements{I<:ID{SimpleID},M<:Element}=Dict{I,M}
"""
    Elements(ms)
    Elements(ms::Pair{I,M}...) where {I<:ID{SimpleID},M<:Element}
    Elements(ms::Element...)

Get the set of elements with similar items merged.
"""
Elements(ms)=Base.dict_with_eltype((K,V)->Dict{K,V},ms,eltype(ms))
function Elements(ms::Pair{I,M}...) where {I<:ID{SimpleID},M<:Element}
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
    show(io::IO,ms::Elements)

Show a set of elements.
"""
function Base.show(io::IO,ms::Elements)
    @printf io "Elements with %s entries:\n" length(ms)
    for m in values(ms)
        @printf io "  %s\n" m
    end
end

"""
    repr(ms::Elements) -> String

Get the repr representation of a set of elements.
"""
function Base.repr(ms::Elements)
    cache=[@sprintf("Elements with %s entries:",length(ms))]
    for m in ms|>values
        push!(cache,@sprintf("  %s",repr(m)))
    end
    return join(cache,"\n")
end

"""
    zero(ms::Elements) -> Nothing
    zero(::Type{<:Elements}) -> Nothing

Get a zero set of elements.

A zero set of elements is defined to be the one with no elements.
"""
Base.zero(ms::Elements)=ms|>typeof|>zero
Base.zero(::Type{<:Elements})=nothing

"""
    ==(ms::Elements,::Nothing) -> Bool
    ==(::Nothing,ms::Elements) -> Bool
    isequal(ms::Elements,::Nothing) -> Bool
    isequal(::Nothing,ms::Elements) -> Bool
"""
Base.:(==)(ms::Elements,::Nothing)=length(ms)==0
Base.:(==)(::Nothing,ms::Elements)=length(ms)==0
Base.isequal(ms::Elements,::Nothing)=length(ms)==0
Base.isequal(::Nothing,ms::Elements)=length(ms)==0

"""
    add!(ms::Elements) -> typeof(ms)
    add!(ms::Elements,::Nothing) -> typeof(ms)
    add!(ms::Elements,m::Number) -> typeof(ms)
    add!(ms::Elements,m::Element) -> typeof(ms)
    add!(ms::Elements,mms::Elements) -> typeof(ms)

Get the inplace addition of elements to a set.
"""
add!(ms::Elements)=ms
add!(ms::Elements,::Nothing)=ms
add!(ms::Elements,m::Number)=add!(ms,convert(scalartype(ms|>valtype),m))
function add!(ms::Elements,m::Element)
    m=convert(ms|>valtype,m)
    old=get(ms,m.id,nothing)
    new=old===nothing ? m : replace(old,value=old.value+m.value)
    abs(new.value)==0.0 ? delete!(ms,m.id) : ms[m.id]=new
    return ms
end
add!(ms::Elements,mms::Elements)=(for m in mms|>values add!(ms,m) end; ms)

"""
    sub!(ms::Elements) -> typeof(ms)
    sub!(ms::Elements,::Nothing) -> typeof(ms)
    add!(ms::Elements,m::Number) -> typeof(ms)
    sub!(ms::Elements,m::Element) -> typeof(ms)
    sub!(ms::Elements,mms::Elements) -> typeof(ms)

Get the inplace subtraction of elements from a set.
"""
sub!(ms::Elements)=ms
sub!(ms::Elements,::Nothing)=ms
sub!(ms::Elements,m::Number)=add!(ms,convert(scalartype(ms|>valtype),-m))
function sub!(ms::Elements,m::Element)
    m=convert(ms|>valtype,m)
    old=get(ms,m.id,nothing)
    new=old==nothing ? -m : replace(old,value=old.value-m.value)
    abs(new.value)==0.0 ? delete!(ms,m.id) : ms[m.id]=new
    return ms
end
sub!(ms::Elements,mms::Elements)=(for m in mms|>values sub!(ms,m) end; ms)

"""
    mul!(ms::Elements,factor::ScalarElement) -> Elements
    mul!(ms::Elements,factor::Number) -> Elements

Get the inplace multiplication of elements with a scalar.
"""
mul!(ms::Elements,factor::ScalarElement)=mul!(ms,factor.value)
function mul!(ms::Elements,factor::Number)
    @assert isa(one(ms|>valtype|>valtype)*factor,ms|>valtype|>valtype) "mul! error: dismatched type, $(ms|>valtype) and $(factor|>typeof)."
    for m in values(ms)
        ms[m.id]=replace(m,value=m.value*factor)
    end
    return ms
end

"""
    div!(ms::Elements,factor::ScalarElement) -> Elements
    div!(ms::Elements,factor::Number) -> Elements

Get the inplace division of element with a scalar.
"""
div!(ms::Elements,factor::ScalarElement)=mul!(ms,1/factor.value)
div!(ms::Elements,factor::Number)=mul!(ms,1/factor)

"""
    +(m::Element) -> typeof(m)
    +(m::Element,::Nothing) -> typeof(m)
    +(::Nothing,m::Element) -> typeof(m)
    +(m::Element,factor::Number) -> Elements
    +(factor::Number,m::Element) -> Elements
    +(m1::Element,m2::Element) -> Elements
    +(ms::Elements) -> typeof(ms)
    +(ms::Elements,::Nothing) -> typeof(ms)
    +(::Nothing,ms::Elements) -> typeof(ms)
    +(ms::Elements,factor::Number) -> Elements
    +(factor::Number,ms::Elements) -> Elements
    +(ms::Elements,m::Element) -> Elements
    +(m::Element,ms::Elements) -> Elements
    +(ms1::Elements,ms2::Elements) -> Elements

Overloaded `+` operator between elements of an algebra over a field.
"""
Base.:+(m::Element)=m
Base.:+(ms::Elements)=ms
Base.:+(m::Element,::Nothing)=m
Base.:+(::Nothing,m::Element)=m
Base.:+(ms::Elements,::Nothing)=ms
Base.:+(::Nothing,ms::Elements)=ms
Base.:+(factor::Number,m::Element)=m+rawtype(m|>typeof)(factor,ID())
Base.:+(m::Element,factor::Number)=m+rawtype(m|>typeof)(factor,ID())
Base.:+(factor::Number,m::ScalarElement)=replace(m,value=m.value+factor)
Base.:+(m::ScalarElement,factor::Number)=replace(m,value=m.value+factor)
Base.:+(m1::ScalarElement,m2::ScalarElement)=replace(m1,value=m1.value+m2.value)
Base.:+(factor::Number,ms::Elements)=ms+rawtype(ms|>valtype)(factor,ID())
Base.:+(ms::Elements,factor::Number)=ms+rawtype(ms|>valtype)(factor,ID())
function Base.:+(m1::Element,m2::Element)
    I=promote_type(m1|>idtype,m2|>idtype)
    M=promote_type(typeof(m1),typeof(m2))
    return add!(Elements{I,M}(m1.id=>m1),m2)
end
Base.:+(m::Element,ms::Elements)=ms+m
function Base.:+(ms::Elements,m::Element)
    I=promote_type(ms|>keytype,m|>idtype)
    M=promote_type(valtype(ms),typeof(m))
    return add!(Elements{I,M}(ms),m)
end
function Base.:+(ms1::Elements,ms2::Elements)
    I=promote_type(ms1|>keytype,ms2|>keytype)
    M=promote_type(valtype(ms1),valtype(ms2))
    return add!(Elements{I,M}(ms1),ms2)
end

"""
    *(m::Element,::Nothing) -> Nothing
    *(::Nothing,m::Element) -> Nothing
    *(factor::Number,m::Element) -> Element
    *(m::Element,factor::Number) -> Element
    *(m1::Element,m2::Element) -> Element
    *(ms::Elements,::Nothing) -> Nothing
    *(::Nothing,ms::Elements) -> Nothing
    *(factor::Number,ms::Elements) -> Elements
    *(ms::Elements,factor::Number) -> Elements
    *(m::Element,ms::Elements) -> Elements
    *(ms::Elements,m::Element) -> Elements
    *(ms1::Elements,ms2::Elements) -> Elements

Overloaded `*` operator for element-scalar multiplications and element-element multiplications of an algebra over a field.
"""
Base.:*(m::Element,::Nothing)=nothing
Base.:*(::Nothing,m::Element)=nothing
Base.:*(ms::Elements,::Nothing)=nothing
Base.:*(::Nothing,ms::Elements)=nothing
Base.:*(factor::ScalarElement,m::Element)=m*factor.value
Base.:*(m::Element,factor::ScalarElement)=m*factor.value
Base.:*(m1::ScalarElement,m2::ScalarElement)=replace(m1,value=m1.value*m2.value)
Base.:*(factor::Number,m::Element)=m*factor
Base.:*(m::Element,factor::Number)=replace(m,value=factor*m.value)
Base.:*(factor::ScalarElement,ms::Elements)=ms*factor.value
Base.:*(ms::Elements,factor::ScalarElement)=ms*factor.value
Base.:*(factor::Number,ms::Elements)=ms*factor
function Base.:*(ms::Elements,factor::Number)
    abs(factor)==0 && return zero(Elements)
    result=Elements{ms|>keytype,mul_promote_type(ms|>valtype,factor|>typeof)}()
    for (id,m) in ms result[id]=m*factor end
    return result
end
Base.:*(m::Element,ms::Elements)=Elements((m*mm for mm in ms|>values)...)
Base.:*(ms::Elements,m::Element)=Elements((mm*m for mm in ms|>values)...)
Base.:*(ms1::Elements,ms2::Elements)=Elements((m1*m2 for m1 in ms1|>values for m2 in ms2|>values)...)
function Base.:*(m1::Element,m2::Element)
    @assert(    m1|>typeof|>nameof==m2|>typeof|>nameof && m1|>typeof|>fieldcount==m2|>typeof|>fieldcount==2,
                "\"*\" error: not implemented between $(m1|>typeof|>nameof) and $(m2|>typeof|>nameof)."
                )
    rawtype(typeof(m1))(m1.value*m2.value,m1.id*m2.id)
end

"""
    -(m::Element) -> typeof(m)
    -(m::Element,::Nothing) -> typeof(m)
    -(::Nothing,m::Element) -> typeof(m)
    -(m::Element,factor::Number) -> Elements
    -(factor::Number,m::Element) -> Elements
    -(m1::Element,m2::Element) -> Elements
    -(ms::Elements) -> typeof(ms)
    -(ms::Elements,::Nothing) -> typeof(ms)
    -(::Nothing,ms::Elements) -> typeof(ms)
    -(ms::Elements,factor::Number) -> Elements
    -(factor::Number,ms::Elements) -> Elements
    -(m::Element,ms::Elements) -> Elements
    -(ms::Elements,m::Element) -> Elements
    -(ms1::Elements,ms2::Elements) -> Elements

Overloaded `-` operator between elements of an algebra over a field.
"""
Base.:-(m::Element)=m*(-1)
Base.:-(ms::Elements)=ms*(-1)
Base.:-(m::Element,::Nothing)=m
Base.:-(::Nothing,m::Element)=-m
Base.:-(ms::Elements,::Nothing)=ms
Base.:-(::Nothing,ms::Elements)=-ms
Base.:-(factor::Number,m::Element)=rawtype(m|>typeof)(factor,ID())-m
Base.:-(m::Element,factor::Number)=m+rawtype(m|>typeof)(-factor,ID())
Base.:-(factor::Number,m::ScalarElement)=replace(m,value=factor-m.value)
Base.:-(m::ScalarElement,factor::Number)=replace(m,value=m.value-factor)
Base.:-(m1::ScalarElement,m2::ScalarElement)=replace(m1,value=m1.value-m2.value)
Base.:-(factor::Number,ms::Elements)=rawtype(ms|>valtype)(factor,ID())-ms
Base.:-(ms::Elements,factor::Number)=ms+rawtype(ms|>valtype)(-factor,ID())
function Base.:-(m1::Element,m2::Element)
    I=promote_type(m1|>idtype,m2|>idtype)
    M=promote_type(typeof(m1),typeof(m2))
    return sub!(Elements{I,M}(m1.id=>m1),m2)
end
function Base.:-(m::Element,ms::Elements)
    I=promote_type(m|>idtype,ms|>keytype)
    M=promote_type(typeof(m),valtype(ms))
    return sub!(Elements{I,M}(m.id=>m),ms)
end
function Base.:-(ms::Elements,m::Element)
    I=promote_type(ms|>keytype,m|>idtype)
    M=promote_type(valtype(ms),typeof(m))
    return sub!(Elements{M|>idtype,M}(ms),m)
end
function Base.:-(ms1::Elements,ms2::Elements)
    I=promote_type(ms1|>keytype,ms2|>keytype)
    M=promote_type(valtype(ms1),valtype(ms2))
    return sub!(Elements{I,M}(ms1),ms2)
end

"""
    /(m::Element,factor::Number) -> Element
    /(m::Element,factor::ScalarElement) -> Element
    /(ms::Elements,factor::Number) -> Elements
    /(ms::Elements,factor::ScalarElement) -> Elements

Overloaded `/` operator for element-scalar division of an algebra over a field.
"""
Base.:/(m::Element,factor::Number)=m*(1/factor)
Base.:/(m::Element,factor::ScalarElement)=m*(1/factor.value)
Base.:/(ms::Elements,factor::Number)=ms*(1/factor)
Base.:/(ms::Elements,factor::ScalarElement)=ms*(1/factor.value)

"""
    ^(m::Element,n::Integer) -> Element
    ^(ms::Elements,n::Integer) -> Elements

Overloaded `^` operator for element-integer power of an algebra over a field.
"""
Base.:^(m::Element,n::Integer)=(@assert n>0 "^ error: non-positive integers are not allowed."; prod(ntuple(i->m,Val(n))))
Base.:^(ms::Elements,n::Integer)=(@assert n>0 "^ error: non-positive integers are not allowed."; prod(ntuple(i->ms,Val(n))))

"""
    ⊗(m::Element,ms::Elements) -> Elements
    ⊗(ms::Elements,m::Element) -> Elements
    ⊗(ms1::Elements,ms2::Elements) -> Elements

Overloaded `⊗` operator for element-element multiplications of an algebra over a field.
"""
⊗(m::Element,ms::Elements)=Elements((m⊗mm for mm in ms|>values)...)
⊗(ms::Elements,m::Element)=Elements((mm⊗m for mm in ms|>values)...)
⊗(ms1::Elements,ms2::Elements)=Elements((m1⊗m2 for m1 in ms1|>values for m2 in ms2|>values)...)

"""
    ⋅(m::Element,ms::Elements) -> Elements
    ⋅(ms::Elements,m::Element) -> Elements
    ⋅(ms1::Elements,ms2::Elements) -> Elements

Overloaded `⋅` operator for element-element multiplications of an algebra over a field.
"""
⋅(m::Element,ms::Elements)=Elements((m⋅mm for mm in ms|>values)...)
⋅(ms::Elements,m::Element)=Elements((mm⋅m for mm in ms|>values)...)
⋅(ms1::Elements,ms2::Elements)=Elements((m1⋅m2 for m1 in ms1|>values for m2 in ms2|>values)...)

"""
    sequence(m::Element,table=nothing) -> NTuple{rank(m),Int}

Get the sequence of the ids of an element according to a table.
"""
@generated sequence(m::Element,table)=Expr(:tuple,[:(get(table,m.id[$i],nothing)) for i=1:rank(m)]...)

"""
    split(m::Element) -> Tuple{Any,Vararg{Element}}

Split an element into the coefficient and a sequence of rank-1 elements.
"""
@generated function Base.split(m::Element)
    @assert m|>fieldnames==(:value,:id) "split error: not supported split of $(nameof(typeof(m)))."
    exprs=[:(m.value)]
    for i=1:rank(m)
        push!(exprs,:(rawtype(typeof(m))(one(m|>valtype),ID(m.id[$i]))))
    end
    return Expr(:tuple,exprs...)
end

"""
    replace(m::Element,pairs::Pair{<:SimpleID,<:Union{Element,Elements}}...) -> Element/Elements
    replace(ms::Elements,pairs::Pair{<:SimpleID,<:Union{Element,Elements}}...) -> Elements

Replace the rank-1 components of an element with new element/elements.
"""
function Base.replace(m::Element,pairs::Pair{<:SimpleID,<:Union{Element,Elements}}...)
    @assert m|>typeof|>fieldnames==(:value,:id) "replace error: not supported replacement of $(nameof(typeof(m)))."
    replacedids=NTuple{length(pairs),idtype(m)|>eltype}(pair.first for pair in pairs)
    ms=split(m)
    result=ms[1]
    for i=1:rank(m)
        index=findfirst(isequal(m.id[i]),replacedids)
        result=result*(isa(index,Int) ? pairs[index].second : ms[i+1])
    end
    return result
end
function Base.replace(ms::Elements,pairs::Pair{<:SimpleID,<:Union{Element,Elements}}...)
    result=elementstype(pairs...)()
    for m in values(ms)
        add!(result,replace(m,pairs...))
    end
    return result
end
@generated function elementstype(ms::Vararg{Pair{<:SimpleID,<:Union{Element,Elements}},N}) where N
    ms=ntuple(i->(fieldtype(ms[i],2)<:Elements) ? fieldtype(ms[i],2)|>valtype : fieldtype(ms[i],2),Val(N))
    @assert mapreduce(m->(m|>fieldnames)==(:value,:id),&,ms) "elementstype error: not supported."
    I=mapreduce(idtype,promote_type,ms)|>eltype
    V=mapreduce(valtype,promote_type,ms)
    M=reduce(typejoin,ms)|>rawtype
    isconcretetype(V) && return Elements{ID{I},M{V,<:ID{I}}}
end

"""
    permute!(result::Elements,m::Element,table=nothing) -> Elements
    permute!(result::Elements,ms::Elements,table=nothing) -> Elements

Permute the ids of an-element/a-set-of-elements to the descending order according to a table, and store the permuted elements in result.
"""
function Base.permute!(result::Elements,m::Element,table=nothing)
    cache=valtype(result)[m]
    while length(cache)>0
        current=pop!(cache)
        pos=elementcommuteposition(sequence(current,table))
        if isa(pos,Nothing)
            add!(result,current)
        else
            left,right=elementleft(current,pos),elementright(current,pos)
            for middle in permute(typeof(current),current.id[pos],current.id[pos+1],table)
                temp=left*middle*right
                temp===nothing || push!(cache,temp)
            end
        end
    end
    return result
end
function Base.permute!(result::Elements,ms::Elements,table=nothing)
    for m in values(ms)
        permute!(result,m,table)
    end
    return result
end
function elementcommuteposition(seqs)
    pos=1
    while pos<length(seqs)
        seqs[pos]<seqs[pos+1] && return pos
        pos+=1
    end
    return nothing
end
elementleft(m::Element,i::Int)=rawtype(typeof(m))(m.value,m.id[1:i-1])
elementright(m::Element,i::Int)=rawtype(typeof(m))(one(m.value),m.id[i+2:end])

"""
    permute(m::Element,table=nothing) -> Elements
    permute(ms::Elements,table=nothing) -> Elements

Permute the ids of an-element/a-set-of-elements to the descending order according to a table.
"""
function permute(m::Element,table=nothing)
    M,V,I=m|>typeof|>rawtype,m|>valtype,m|>idtype|>eltype
    permute!(Elements{ID{I},M{V,<:ID{I}}}(),m,table)
end
function permute(ms::Elements,table=nothing)
    M,V,I=ms|>valtype|>rawtype,ms|>valtype|>valtype,ms|>keytype|>eltype
    permute!(Elements{ID{I},M{V,<:ID{I}}}(),ms,table)
end

end #module
