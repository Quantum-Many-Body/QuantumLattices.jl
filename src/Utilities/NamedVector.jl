module NamedVector

using Printf: @printf
using ..Factory: Inference,TypeFactory,FunctionFactory,addargs!,extendbody!
using ..Factory: MixEscaped,Escaped,UnEscaped
using ..Utilities: efficientoperations

export AbstractNamedVector,@namedvector
export HomoNamedVector,@homonamedvector

"""
    AbstractNamedVector

Abstract type for all named vectors.
"""
abstract type AbstractNamedVector end

"""
    getindex(nv::AbstractNamedVector,index::Int)

Get the value by the `[]` syntax.
"""
Base.getindex(nv::AbstractNamedVector,index::Int)=getfield(nv,index)

"""
    setindex!(nv::AbstractNamedVector,value,index::Int)

Set the value by the `[]` syntax if mutable.
"""
Base.setindex!(nv::AbstractNamedVector,value,index::Int)=setfield!(nv,index,value)

"""
    ==(nv1::AbstractNamedVector,nv2::AbstractNamedVector) -> Bool
    isequal(nv1::AbstractNamedVector,nv2::AbstractNamedVector) -> Bool

Overloaded equivalent operator. Two named vector are equal to each other if and only if their keys as well as their values are equal to each other.
!!! note
    It is not necessary for two named vectors to be of the same concrete type to be equal to each other.
"""
Base.:(==)(nv1::AbstractNamedVector,nv2::AbstractNamedVector)=nv1|>keys==nv2|>keys && nv1|>values==nv2|>values
Base.isequal(nv1::AbstractNamedVector,nv2::AbstractNamedVector)=isequal(nv1|>keys,nv2|>keys) && isequal(nv1|>values,nv2|>values)

"""
    <(nv1::AbstractNamedVector,nv2::AbstractNamedVector) -> Bool
    isless(nv1::AbstractNamedVector,nv2::AbstractNamedVector) -> Bool

Compare two named vectors and judge whether the first is less than the second.
"""
Base.:<(nv1::AbstractNamedVector,nv2::AbstractNamedVector) = <(efficientoperations,nv1,nv2)
Base.isless(nv1::AbstractNamedVector,nv2::AbstractNamedVector)=isless(efficientoperations,nv1,nv2)

"""
    show(io::IO,nv::AbstractNamedVector)

Show a concrete `AbstractNamedVector`.
"""
Base.show(io::IO,nv::AbstractNamedVector)=@printf io "%s(%s)" nv|>typeof|>nameof join(repr.(nv|>values),',')

"""
    hash(nv::AbstractNamedVector,h::UInt)

Hash a concrete `AbstractNamedVector`.
"""
Base.hash(nv::AbstractNamedVector,h::UInt)=hash(nv|>values,h)

"""
    convert(::Type{Tuple},nv::AbstractNamedVector) -> Tuple
    convert(::Type{NV},nv::Tuple) where NV<:AbstractNamedVector -> NV

Convert a named vector to tuple and vice versa.
"""
@generated Base.convert(::Type{Tuple},nv::AbstractNamedVector)=Expr(:tuple,(:(getfield(nv,$i)) for i=1:(nv|>fieldcount))...)
function Base.convert(::Type{NV},nv::Tuple) where NV<:AbstractNamedVector
    @assert NV|>fieldcount==nv|>length "convert error: dismatched length between $NV($(NV|>fieldcount)) and input tuple($(nv|>length))."
    return NV(nv...)
end

"""
    length(::Type{NV}) where NV<:AbstractNamedVector -> Int
    length(nv::AbstractNamedVector) -> Int

Get the length of a concrete `AbstractNamedVector`.
"""
Base.length(::Type{NV}) where NV<:AbstractNamedVector=NV|>fieldcount
Base.length(nv::AbstractNamedVector)=nv|>typeof|>length

"""
    zero(::Type{NV}) where NV<:AbstractNamedVector -> NV
    zero(nv::AbstractNamedVector) -> typeof(nv)

Get a concrete `AbstractNamedVector` with all values being zero.
"""
@generated Base.zero(::Type{NV}) where NV<:AbstractNamedVector=(zeros=(zero(fieldtype(NV,i)) for i=1:(NV|>fieldcount));:(NV($(zeros...))))
Base.zero(nv::AbstractNamedVector)=nv|>typeof|>zero

"""
    iterate(nv::AbstractNamedVector,state=1)
    iterate(rv::Iterators.Reverse{<:AbstractNamedVector},state=length(rv.itr))

Iterate or reversely iterate over the values of a concrete `AbstractNamedVector`.
"""
Base.iterate(nv::AbstractNamedVector,state=1)=state>length(nv) ? nothing : (getfield(nv,state),state+1)
Base.iterate(rv::Iterators.Reverse{<:AbstractNamedVector},state=length(rv.itr))=state<1 ? nothing : (getfield(rv.itr,state),state-1)

"""
    keys(nv::AbstractNamedVector) -> NTuple(nv|>fieldcount,Symbol)

Iterate over the names.
"""
Base.keys(nv::AbstractNamedVector)=nv|>typeof|>fieldnames

"""
    values(nv::AbstractNamedVector) -> Tuple

Iterate over the values.
"""
Base.values(nv::AbstractNamedVector)=convert(Tuple,nv)

"""
    pairs(nv::AbstractNamedVector)

Iterate over the `name=>value` pairs.
"""
Base.pairs(nv::AbstractNamedVector)=Base.Generator(=>,keys(nv),values(nv))

"""
    replace(nv::AbstractNamedVector;kwargs...) -> typeof(nv)

Return a copy of a concrete `AbstractNamedVector` with some of the field values replaced by the keyword arguments.
"""
Base.replace(nv::AbstractNamedVector;kwargs...)=replace(efficientoperations,nv;kwargs...)

"""
    map(f,nvs::NV...) where NV<:AbstractNamedVector -> NV

Apply function `f` elementwise on the input named vectors.
"""
@generated function Base.map(f,nvs::NV...) where NV<:AbstractNamedVector
    exprs=Vector{Expr}(undef,NV|>fieldcount)
    for i=1:length(exprs)
        tmp=[:(nvs[$j][$i]) for j=1:length(nvs)]
        exprs[i]=:(f($(tmp...)))
    end
    return :(($NV)($(exprs...)))
end

"""
    @namedvector structdef::Expr

Decorate a "raw" struct to be a subtype of `AbstractNamedVector`. Here, "raw" means that the input struct has no explicit supertype and no inner constructors.
"""
macro namedvector(structdef::Expr)
    tf=TypeFactory(structdef)
    @assert tf.supertype==Inference(:Any) "@namedvector error: no explicit supertype except `Any` is allowed."
    @assert length(tf.constructors)==0 "@namedvector error: no inner constructor is allowed."
    tf.supertype=Inference(:AbstractNamedVector)
    paramnames=tuple((param.name for param in tf.params)...)
    fieldnames=tuple((field.name for field in tf.fields)...)
    fldnm=FunctionFactory(name=:(Base.fieldnames))
    addargs!(fldnm,:(::Type{<:$(tf.name)}))
    extendbody!(fldnm,Expr(:tuple,QuoteNode.(fieldnames)...))
    structdef=tf(MixEscaped(Escaped(tf.name),UnEscaped(paramnames...,:AbstractNamedVector)))
    fldnmdef=fldnm(MixEscaped(Escaped(:tuple)))
    return Expr(:block,:(Base.@__doc__($structdef)),fldnmdef)
end

"""
    HomoNamedVector{T}

Abstract type for all homogeneous named vectors.
"""
abstract type HomoNamedVector{T} <: AbstractNamedVector end

"""
    eltype(::Type{NV}) where NV<:HomoNamedVector{T} where T
    eltype(nv::HomoNamedVector)

Get the type parameter of a concrete `HomoNamedVector`.
"""
Base.eltype(::Type{<:HomoNamedVector{T}}) where T=T
Base.eltype(nv::HomoNamedVector)=nv|>typeof|>eltype

"""
    @homonamedvector typename fieldnames dtype::Union{Expr,Symbol}=:nothing mutable::Union{Expr,Bool}=false

Construct a concrete homogeneous named vector with the type name being `typename` and the fieldnames specified by `fieldnames`, and optionally, the type parameters specified by `dtype`.`mutable` can be used as a keyword argument to determine whether the concrete type is mutable.
"""
macro homonamedvector(typename,fieldnames,dtype::Union{Expr,Symbol}=:nothing,mutable::Union{Expr,Bool}=false)
    typename=Symbol(typename)
    fieldnames=tuple(eval(fieldnames)...)
    @assert all(isa(name,Symbol) for name in fieldnames) "homonamedvector error: every field name should be a `Symbol`."
    isa(dtype,Expr) && (@assert (dtype.head==:(<:) && dtype.args|>length==1) "homonamedvector error: wrong `dtype`.")
    dname,dscope=isa(dtype,Expr) ? (:T,dtype.args[1]) : dtype==:nothing ? (:T,:Any) : (dtype,:concrete)
    if isa(mutable,Expr)
        @assert mutable.head==:(=) && mutable.args[1]==:mutable && isa(mutable.args[2],Bool) "homonamedvector error: wrong `mutable`."
        mutable=mutable.args[2]
    end
    if dscope==:concrete
        new=:($(esc(typename)))
        super=:(HomoNamedVector{$(esc(dname))})
        body=(:($field::$(esc(dname))) for field in fieldnames)
    else
        new=:($(esc(typename)){$dname<:$(esc(dscope))})
        super=:(HomoNamedVector{$dname})
        body=(:($field::$dname) for field in fieldnames)
    end
    structdef=Expr(:struct,mutable,Expr(:<:,new,super),Expr(:block,body...))
    functions=:(Base.fieldnames(::Type{<:$(esc(typename))})=$fieldnames)
    return Expr(:block,:(Base.@__doc__($structdef)),functions)
end

end #module
