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

abstract type SimpleID <: AbstractNamedVector end
Base.eltype(id::SimpleID)=id|>typeof|>eltype
Base.eltype(::Type{I}) where I<:SimpleID=I
Base.length(id::SimpleID)=id|>typeof|>length
Base.length(::Type{<:SimpleID})=1
Base.iterate(id::SimpleID,state::Integer=1)=state>1 ? nothing : (id,state+1)
Base.iterate(rv::Iterators.Reverse{<:SimpleID},state::Integer=1)=state<1 ? nothing : (rv.itr,state-1)
rank(::Type{<:SimpleID})=1
rank(id::SimpleID)=id|>typeof|>rank

struct CompositeID{N,I<:SimpleID} <: CompositeNTuple{N,I}
    contents::NTuple{N,I}
    CompositeID(ids::NTuple{N,SimpleID}) where N=new{N,ids|>eltype}(ids)
end
CompositeID(ids::SimpleID...)=CompositeID(ids)
Base.show(io::IO,cid::CompositeID)=@printf io "%s(%s)" cid|>typeof|>nameof join(cid,",")
Base.typejoin(SID::Type{<:SimpleID},CID::Type{<:CompositeID})=Union{SID,CID}
Base.typejoin(CID::Type{<:CompositeID},SID::Type{<:SimpleID})=Union{CID,SID}
rank(::Type{<:CompositeID{N,I}}) where {N,I}=N
rank(id::CompositeID)=id|>typeof|>rank

const ID=Union{<:SimpleID,<:CompositeID}
⊗(sid1::SimpleID,sid2::SimpleID)=CompositeID(sid1,sid2)
⊗(sid::SimpleID,cid::CompositeID)=CompositeID(sid,convert(Tuple,cid)...)
⊗(cid::CompositeID,sid::SimpleID)=CompositeID(convert(Tuple,cid)...,sid)
⊗(cid1::CompositeID,cid2::CompositeID)=CompositeID(convert(Tuple,cid1)...,convert(Tuple,cid2)...)
function Base.isless(cid1::ID,cid2::ID)
    r1,r2=cid1|>rank,cid2|>rank
    r1<r2 ? true : r1>r2 ? false : isless(convert(Tuple,cid1),convert(Tuple,cid2))
end
function Base.:<(cid1::ID,cid2::ID)
    r1,r2=cid1|>rank,cid2|>rank
    r1<r2 ? true : r1>r2 ? false : convert(Tuple,cid1)<convert(Tuple,cid2)
end


const SimpleVectorSpace{I<:ID}=Vector{I}
SimpleVectorSpace(ids::ID...)=collect(ids)
const CompositeVectorSpace{SVS<:SimpleVectorSpace}=Vector{SVS}
CompositeVectorSpace(svses::SimpleVectorSpace...)=collect(svses)
const VectorSpace=Union{SimpleVectorSpace,CompositeVectorSpace}
VectorSpace(ids::ID...)=SimpleVectorSpace(ids...)
VectorSpace(svses::SimpleVectorSpace...)=CompositeVectorSpace(svses...)

dimension(svs::SimpleVectorSpace)=length(svs)
dimension(cvs::CompositeVectorSpace)=prod(length(svs) for svs in cvs)

⊕(id1::I,id2::I) where {I<:ID}=push!(SimpleVectorSpace{I}(),id1,id2)
⊕(id::I,svs::SimpleVectorSpace{I}) where {I<:ID}=append!(push!(SimpleVectorSpace{I}(),id),svs)
⊕(svs::SimpleVectorSpace{I},id::I) where {I<:ID}=push!(append!(SimpleVectorSpace{I}(),svs),id)
⊕(svs1::SVS,svs2::SVS) where {SVS<:SimpleVectorSpace}=append!(append!(SimpleVectorSpace{SVS|>eltype}(),svs1),svs2)

⊗(svs1::SimpleVectorSpace,svs2::SimpleVectorSpace)=push!(CompositeVectorSpace{typejoin(svs1|>typeof,svs2|>typeof)}(),svs1,svs2)
⊗(svs::SimpleVectorSpace,cvs::CompositeVectorSpace)=append!(push!(CompositeVectorSpace{typejoin(svs|>typeof,cvs|>typeof|>eltype)}(),svs),cvs)
⊗(cvs::CompositeVectorSpace,svs::SimpleVectorSpace)=push!(append!(CompositeVectorSpace{typejoin(cvs|>typeof|>eltype,svs|>typeof)}(),cvs),svs)
⊗(cvs1::CompositeVectorSpace,cvs2::CompositeVectorSpace)=append!(append!(CompositeVectorSpace{typejoin(cvs1|>typeof|>eltype,cvs2|>typeof|>eltype)}(),cvs1),cvs2)


abstract type Element{V<:Number,I<:ID} end
Base.valtype(::Type{<:Element{V,I}}) where {V,I}=V
Base.valtype(m::Element)=m|>typeof|>valtype
idtype(::Type{<:Element{V,I}}) where {V,I}=I
idtype(m::Element)=m|>typeof|>idtype
rank(::Type{<:Element{V,I}}) where {V,I}=I|>rank
rank(m::Element)=m|>typeof|>rank
id(m::Element)=m.id

const Elements{I<:ID,M<:Element}=Dict{I,M}
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


Base.zero(ms::Elements)=ms|>typeof|>zero
Base.zero(::Type{Elements{I,M}}) where {I,M}=Elements{I,M}()

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

sub!(ms::Elements)=ms
function sub!(ms::Elements,m::Element)
    @assert ms|>valtype==m|>typeof "sub! error: dismatched type, $(ms|>valtype) and $(m|>typeof)."
    mid=m|>id
    ms[mid]=haskey(ms,mid) ? (m|>typeof)(ms[mid].value-m.value,(getfield(m,i) for i=2:(m|>typeof|>fieldcount))...) : -m
    abs(ms[mid].value)==0.0 && delete!(ms,mid)
    ms
end
sub!(ms::Elements,mms::Elements)=(for m in mms|>values sub!(ms,m) end; ms)

Base.:+(m::Element)=m
Base.:+(ms::Elements)=ms
Base.:+(ms::Elements,m::Element)=m+ms
Base.:+(m1::Element,m2::Element)=add!(Elements{typejoin(m1|>idtype,m2|>idtype),typejoin(m1|>typeof,m2|>typeof)}(m1|>id=>m1),m2)
Base.:+(m::Element,ms::Elements)=add!(Elements{typejoin(m|>idtype,ms|>keytype),typejoin(m|>typeof,ms|>valtype)}(ms),m)
Base.:+(ms1::Elements,ms2::Elements)=add!(Elements{typejoin(ms1|>keytype,ms2|>keytype),typejoin(ms1|>valtype,ms2|>valtype)}(ms1),ms2)

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

Base.:-(m::Element)=m*(-1)
Base.:-(ms::Elements)=ms*(-1)
Base.:-(m1::Element,m2::Element)=sub!(Elements{typejoin(m1|>idtype,m2|>idtype),typejoin(m1|>typeof,m2|>typeof)}(m1|>id=>m1),m2)
Base.:-(m::Element,ms::Elements)=sub!(Elements{typejoin(m|>idtype,ms|>keytype),typejoin(m|>typeof,ms|>valtype)}(m|>id=>m),ms)
Base.:-(ms::Elements,m::Element)=sub!(Elements{typejoin(m|>idtype,ms|>keytype),typejoin(m|>typeof,ms|>valtype)}(ms),m)
Base.:-(ms1::Elements,ms2::Elements)=sub!(Elements{typejoin(ms1|>keytype,ms2|>keytype),typejoin(ms1|>valtype,ms2|>valtype)}(ms1),ms2)

Base.:/(m::Element,factor::Number)=m*(1/factor)
Base.:/(ms::Elements,factor::Number)=ms*(1/factor)

end #module
