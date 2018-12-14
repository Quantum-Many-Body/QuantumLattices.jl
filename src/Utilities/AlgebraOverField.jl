module AlgebraOverField

using Printf: @printf
using ..NamedVector: AbstractNamedVector
using ..CompositeStructure: CompositeNTuple

import ..Utilities: rank,⊕,⊗,dimension

export SimpleID,CompositeID,ID
export SimpleVectorSpace,CompositeVectorSpace,VectorSpace
export rank,⊕,⊗,dimension
export Element,Elements
export idtype,id,add!,sub!

"""
    SimpleID <: AbstractNamedVector

A simple id is the id of a single basis of a vector space or algebra over a field.
"""
abstract type SimpleID <: AbstractNamedVector end

"""
    eltype(id::SimpleID)
    eltype(::Type{I}) where I

Get the eltype of a simple id, which is defined to be the type of itself.
"""
Base.eltype(id::SimpleID)=id|>typeof|>eltype
Base.eltype(::Type{I}) where I<:SimpleID=I

"""
    length(id::SimpleID) -> Int
    length(::Type{<:SimpleID}) -> Int

Get the length of a simple id, which is defined to be 1.
"""
Base.length(id::SimpleID)=id|>typeof|>length
Base.length(::Type{<:SimpleID})=1

"""
    iterate(id::SimpleID,state::Integer=1)
    iterate(rv::Iterators.Reverse{<:SimpleID},state::Integer=1)

Iterate over a simple id.

The iteration is defined to give itself.
"""
Base.iterate(id::SimpleID,state::Integer=1)=state>1 ? nothing : (id,state+1)
Base.iterate(rv::Iterators.Reverse{<:SimpleID},state::Integer=1)=state<1 ? nothing : (rv.itr,state-1)

"""
    rank(::Type{<:SimpleID}) -> Int
    rank(id::SimpleID) -> Int

Get the rank of a simple id, which is defined to be 1.
"""
rank(::Type{<:SimpleID})=1
rank(id::SimpleID)=id|>typeof|>rank

"""
    CompositeID(ids::NTuple{N,SimpleID}) where N
    CompositeID(ids::SimpleID...)

A composite id is the id of the multiplication of bases of an algebra over a field.
"""
struct CompositeID{N,I<:SimpleID} <: CompositeNTuple{N,I}
    contents::NTuple{N,I}
    CompositeID(ids::NTuple{N,SimpleID}) where N=new{N,ids|>eltype}(ids)
end
CompositeID(ids::SimpleID...)=CompositeID(ids)

"""
    show(io::IO,cid::CompositeID)

Show a composite id.
"""
Base.show(io::IO,cid::CompositeID)=@printf io "%s(%s)" cid|>typeof|>nameof join(cid,",")

"""
    hash(cid::CompositeID,h::UInt)

Hash a composite id.
"""
Base.hash(cid::CompositeID,h::UInt)=hash(convert(Tuple,cid),h)

"""
    rank(::Type{<:CompositeID{N,I}}) where {N,I} -> Int
    rank(id::CompositeID) -> Int

Get the rank of a composite id.
"""
rank(::Type{<:CompositeID{N,I}}) where {N,I}=N
rank(id::CompositeID)=id|>typeof|>rank

"""
    typejoin(SID::Type{<:SimpleID},CID::Type{<:CompositeID})
    typejoin(CID::Type{<:CompositeID},SID::Type{<:SimpleID})

Get the type join of a simple id and a composite id.
"""
Base.typejoin(SID::Type{<:SimpleID},CID::Type{<:CompositeID})=Union{SID,CID}
Base.typejoin(CID::Type{<:CompositeID},SID::Type{<:SimpleID})=Union{CID,SID}

"""
    ID

The id system of an algebra over a field.

It is a type alias for `Union{<:SimpleID,<:CompositeID}`.
"""
const ID=Union{<:SimpleID,<:CompositeID}

"""
    ⊗(sid1::SimpleID,sid2::SimpleID) -> CompositeID
    ⊗(sid::SimpleID,cid::CompositeID) -> CompositeID
    ⊗(cid::CompositeID,sid::SimpleID) -> CompositeID
    ⊗(cid1::CompositeID,cid2::CompositeID) -> CompositeID

Get the direct product of the id system.
"""
⊗(sid1::SimpleID,sid2::SimpleID)=CompositeID(sid1,sid2)
⊗(sid::SimpleID,cid::CompositeID)=CompositeID(sid,convert(Tuple,cid)...)
⊗(cid::CompositeID,sid::SimpleID)=CompositeID(convert(Tuple,cid)...,sid)
⊗(cid1::CompositeID,cid2::CompositeID)=CompositeID(convert(Tuple,cid1)...,convert(Tuple,cid2)...)

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
    SimpleVectorSpace(ids::ID...)

The vector space spanned by a set of bases specified by their ids.
"""
struct SimpleVectorSpace{I<:ID,N} <: CompositeNTuple{N,I}
    contents::NTuple{N,I}
end
SimpleVectorSpace(ids::ID...)=SimpleVectorSpace(ids)

"""
    CompositeVectorSpace(svses::SimpleVectorSpace...)

The vector space spanned by the direct product of simple vector spaces.
"""
struct CompositeVectorSpace{SVS<:SimpleVectorSpace,N} <: CompositeNTuple{N,SVS}
    contents::NTuple{N,SVS}
end
CompositeVectorSpace(svses::SimpleVectorSpace...)=CompositeVectorSpace(svses)

"""
    VectorSpace

The corresponding vector space of an algebra over a field.

Alias for `Union{SimpleVectorSpace,CompositeVectorSpace}`.
"""
const VectorSpace=Union{SimpleVectorSpace,CompositeVectorSpace}
"""
    VectorSpace(ids::ID...) -> SimpleVectorSpace
    VectorSpace(svses::SimpleVectorSpace...) -> CompositeVectorSpace

Get the corresponding vector space of an algebra over a field.
"""
VectorSpace(ids::ID...)=SimpleVectorSpace(ids...)
VectorSpace(svses::SimpleVectorSpace...)=CompositeVectorSpace(svses...)

"""
    dimension(svs::SimpleVectorSpace) -> Int
    dimension(cvs::CompositeVectorSpace) -> Int

Get the dimension of a vector space.
"""
dimension(svs::SimpleVectorSpace)=length(svs)
dimension(cvs::CompositeVectorSpace)=prod(length(svs) for svs in cvs)

"""
    ⊕(id1::I,id2::I) where {I<:ID} -> SimpleVectorSpace{I}
    ⊕(id::I,svs::SimpleVectorSpace{I}) where {I<:ID} -> SimpleVectorSpace{I}
    ⊕(svs::SimpleVectorSpace{I},id::I) where {I<:ID} -> SimpleVectorSpace{I}
    ⊕(svs1::SVS,svs2::SVS) where {SVS<:SimpleVectorSpace} -> SVS

Get the direct sum of bases or simple vector spaces.
"""
⊕(id1::I,id2::I) where {I<:ID}=SimpleVectorSpace(id1,id2)
⊕(id::I,svs::SimpleVectorSpace{I}) where {I<:ID}=SimpleVectorSpace(id,convert(Tuple,svs)...)
⊕(svs::SimpleVectorSpace{I},id::I) where {I<:ID}=SimpleVectorSpace(convert(Tuple,svs)...,id)
⊕(svs1::SVS,svs2::SVS) where {SVS<:SimpleVectorSpace}=SimpleVectorSpace(convert(Tuple,svs1)...,convert(Tuple,svs2)...)

"""
    ⊗(svs1::SimpleVectorSpace,svs2::SimpleVectorSpace) -> CompositeVectorSpace
    ⊗(svs::SimpleVectorSpace,cvs::CompositeVectorSpace) -> CompositeVectorSpace
    ⊗(cvs::CompositeVectorSpace,svs::SimpleVectorSpace) -> CompositeVectorSpace
    ⊗(cvs1::CompositeVectorSpace,cvs2::CompositeVectorSpace) -> CompositeVectorSpace

Get the direct product of simple vector spaces or composite vector spaces.
"""
⊗(svs1::SimpleVectorSpace,svs2::SimpleVectorSpace)=CompositeVectorSpace(svs1,svs2)
⊗(svs::SimpleVectorSpace,cvs::CompositeVectorSpace)=CompositeVectorSpace(svs,convert(Tuple,cvs)...)
⊗(cvs::CompositeVectorSpace,svs::SimpleVectorSpace)=CompositeVectorSpace(convert(Tuple,cvs)...,svs)
⊗(cvs1::CompositeVectorSpace,cvs2::CompositeVectorSpace)=CompositeVectorSpace(convert(Tuple,cvs1)...,convert(Tuple,cvs2)...)

"""
    Element{V<:Number,I<:ID}

An element of an algebra over a field.
"""
abstract type Element{V<:Number,I<:ID} end

"""
    valtype(::Type{<:Element{V,I}}) where {V,I}
    valtype(m::Element)

Get the type of the value of an element.

The result is also the type of the field over which the algebra is defined.
"""
Base.valtype(::Type{<:Element{V,I}}) where {V,I}=V
Base.valtype(m::Element)=m|>typeof|>valtype

"""
    idtype(::Type{<:Element{V,I}}) where {V,I}
    idtype(m::Element)

The type of the id of an element.
"""
idtype(::Type{<:Element{V,I}}) where {V,I}=I
idtype(m::Element)=m|>typeof|>idtype

"""
    rank(::Type{<:Element{V,I}}) where {V,I} -> Int
    rank(m::Element) -> Int

Get the rank of an element.
"""
rank(::Type{<:Element{V,I}}) where {V,I}=I|>rank
rank(m::Element)=m|>typeof|>rank

"""
    id(m::Element) -> ID

Get the id of an element.
"""
id(m::Element)=m.id


"""
    Elements{I<:ID,M<:Element} -> AbstractDict{I,M}

An set of elements of an algebra over a field.

Similar iterms are automatically merged thanks to the id system.
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
    add!(ms::Elements,mms::Element...) -> typeof(ms)
    add!(ms::Elements,mses::Elements...) -> typeof(ms)

Get the inplace addition of elements to a set.
"""
add!(ms::Elements)=ms
function add!(ms::Elements,mms::Element...)
    mms|>length>0 && @assert ms|>valtype==mms|>eltype "add! error: dismatched type, $(ms|>valtype) and $(mms|>eltype)."
    for mm in mms
        mid=mm|>id
        ms[mid]=haskey(ms,mid) ? (mm|>typeof)(ms[mid].value+mm.value,(getfield(mm,i) for i=2:(mm|>typeof|>fieldcount))...) : mm
        abs(ms[mid].value)==0.0 && delete!(ms,mid)
    end
    ms
end
add!(ms::Elements,mses::Elements...)=(for mms in mses for m in mms|>values add!(ms,m) end end; ms)

"""
    sub!(ms::Elements) -> typeof(ms) -> typeof(ms)
    sub!(ms::Elements,m::Element) -> typeof(ms)
    sub!(ms::Elements,mms::Elements) -> typeof(ms)

Get the inplace subtraction of elements from a set.
"""
sub!(ms::Elements)=ms
function sub!(ms::Elements,m::Element)
    @assert ms|>valtype==m|>typeof "sub! error: dismatched type, $(ms|>valtype) and $(m|>typeof)."
    mid=m|>id
    ms[mid]=haskey(ms,mid) ? (m|>typeof)(ms[mid].value-m.value,(getfield(m,i) for i=2:(m|>typeof|>fieldcount))...) : -m
    abs(ms[mid].value)==0.0 && delete!(ms,mid)
    ms
end
sub!(ms::Elements,mms::Elements)=(for m in mms|>values sub!(ms,m) end; ms)

"""
    +(m::Element)
    +(ms::Elements)
    +(ms::Elements,m::Element)
    +(m1::Element,m2::Element)
    +(m::Element,ms::Elements)
    +(ms1::Elements,ms2::Elements)

Overloaded `+` operator between elements of an algebra over a field.
"""
Base.:+(m::Element)=m
Base.:+(ms::Elements)=ms
Base.:+(ms::Elements,m::Element)=m+ms
Base.:+(m1::Element,m2::Element)=add!(Elements{typejoin(m1|>idtype,m2|>idtype),typejoin(m1|>typeof,m2|>typeof)}(m1|>id=>m1),m2)
Base.:+(m::Element,ms::Elements)=add!(Elements{typejoin(m|>idtype,ms|>keytype),typejoin(m|>typeof,ms|>valtype)}(ms),m)
Base.:+(ms1::Elements,ms2::Elements)=add!(Elements{typejoin(ms1|>keytype,ms2|>keytype),typejoin(ms1|>valtype,ms2|>valtype)}(ms1),ms2)

"""
    *(factor::Number,m::Element)
    *(factor::Number,ms::Elements)
    *(m::Element,factor::Number)
    *(ms::Elements,factor::Number)
    *(m::Element,ms::Elements)
    *(ms::Elements,m::Element)
    *(ms1::Elements,ms2::Elements)
    *(m1::Element,m2::Element)

Overloaded `*` operator for element-scalar multiplications and element-element multiplications of an algebra over a field.
"""
Base.:*(factor::Number,m::Element)=m*factor
Base.:*(factor::Number,ms::Elements)=ms*factor
Base.:*(m::Element,factor::Number)=typeof(m).name.wrapper(m.value*factor,(getfield(m,i) for i=2:(m|>typeof|>fieldcount))...)
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
    -(m::Element)
    -(ms::Elements)
    -(m1::Element,m2::Element)
    -(m::Element,ms::Elements)
    -(ms::Elements,m::Element)
    -(ms1::Elements,ms2::Elements)

Overloaded `-` operator between elements of an algebra over a field.
"""
Base.:-(m::Element)=m*(-1)
Base.:-(ms::Elements)=ms*(-1)
Base.:-(m1::Element,m2::Element)=sub!(Elements{typejoin(m1|>idtype,m2|>idtype),typejoin(m1|>typeof,m2|>typeof)}(m1|>id=>m1),m2)
Base.:-(m::Element,ms::Elements)=sub!(Elements{typejoin(m|>idtype,ms|>keytype),typejoin(m|>typeof,ms|>valtype)}(m|>id=>m),ms)
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
