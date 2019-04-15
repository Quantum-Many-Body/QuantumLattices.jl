module NamedVectors

using Printf: @printf
using ..Factories: Inference,TypeFactory,FunctionFactory,addargs!,extendbody!
using ..Factories: MixEscaped,Escaped,UnEscaped
using ..TypeTraits: efficientoperations

export NamedVector,@namedvector
export HomoNamedVector,@homonamedvector

"""
    NamedVector

Abstract type for all named vectors.
"""
abstract type NamedVector end

"""
    getindex(nv::NamedVector,index::Int)

Get the value by the `[]` syntax.
"""
Base.getindex(nv::NamedVector,index::Int)=getfield(nv,index)

"""
    setindex!(nv::NamedVector,value,index::Int)

Set the value by the `[]` syntax if mutable.
"""
Base.setindex!(nv::NamedVector,value,index::Int)=setfield!(nv,index,value)

"""
    ==(nv1::NamedVector,nv2::NamedVector) -> Bool
    isequal(nv1::NamedVector,nv2::NamedVector) -> Bool

Overloaded equivalent operator. Two named vector are equal to each other if and only if their keys as well as their values are equal to each other.
!!! note
    It is not necessary for two named vectors to be of the same concrete type to be equal to each other.
"""
Base.:(==)(nv1::NamedVector,nv2::NamedVector)=nv1|>keys==nv2|>keys && nv1|>values==nv2|>values
Base.isequal(nv1::NamedVector,nv2::NamedVector)=isequal(nv1|>keys,nv2|>keys) && isequal(nv1|>values,nv2|>values)

"""
    <(nv1::NamedVector,nv2::NamedVector) -> Bool
    isless(nv1::NamedVector,nv2::NamedVector) -> Bool

Compare two named vectors and judge whether the first is less than the second.
"""
Base.:<(nv1::NamedVector,nv2::NamedVector) = <(efficientoperations,nv1,nv2)
Base.isless(nv1::NamedVector,nv2::NamedVector)=isless(efficientoperations,nv1,nv2)

"""
    show(io::IO,nv::NamedVector)

Show a concrete `NamedVector`.
"""
Base.show(io::IO,nv::NamedVector)=@printf io "%s(%s)" nv|>typeof|>nameof join(repr.(nv|>values),',')

"""
    hash(nv::NamedVector,h::UInt)

Hash a concrete `NamedVector`.
"""
Base.hash(nv::NamedVector,h::UInt)=hash(nv|>values,h)

"""
    convert(::Type{Tuple},nv::NamedVector) -> Tuple
    convert(::Type{NV},nv::Tuple) where NV<:NamedVector -> NV

Convert a named vector to tuple and vice versa.
"""
@generated Base.convert(::Type{Tuple},nv::NamedVector)=Expr(:tuple,(:(getfield(nv,$i)) for i=1:(nv|>fieldcount))...)
function Base.convert(::Type{NV},nv::Tuple) where NV<:NamedVector
    @assert NV|>fieldcount==nv|>length "convert error: dismatched length between $NV($(NV|>fieldcount)) and input tuple($(nv|>length))."
    return NV(nv...)
end

"""
    length(::Type{NV}) where NV<:NamedVector -> Int
    length(nv::NamedVector) -> Int

Get the length of a concrete `NamedVector`.
"""
Base.length(::Type{NV}) where NV<:NamedVector=NV|>fieldcount
Base.length(nv::NamedVector)=nv|>typeof|>length

"""
    zero(::Type{NV}) where NV<:NamedVector -> NV
    zero(nv::NamedVector) -> typeof(nv)

Get a concrete `NamedVector` with all values being zero.
"""
@generated Base.zero(::Type{NV}) where NV<:NamedVector=(zeros=(zero(fieldtype(NV,i)) for i=1:(NV|>fieldcount));:(NV($(zeros...))))
Base.zero(nv::NamedVector)=nv|>typeof|>zero

"""
    iterate(nv::NamedVector,state=1)
    iterate(rv::Iterators.Reverse{<:NamedVector},state=length(rv.itr))

Iterate or reversely iterate over the values of a concrete `NamedVector`.
"""
Base.iterate(nv::NamedVector,state=1)=state>length(nv) ? nothing : (getfield(nv,state),state+1)
Base.iterate(rv::Iterators.Reverse{<:NamedVector},state=length(rv.itr))=state<1 ? nothing : (getfield(rv.itr,state),state-1)

"""
    keys(nv::NamedVector) -> NTuple(nv|>fieldcount,Symbol)
    values(nv::NamedVector) -> Tuple
    pairs(nv::NamedVector) -> Base.Generator

Iterate over the names.
"""
Base.keys(nv::NamedVector)=nv|>typeof|>fieldnames
Base.values(nv::NamedVector)=convert(Tuple,nv)
Base.pairs(nv::NamedVector)=Base.Generator(=>,keys(nv),values(nv))

"""
    replace(nv::NamedVector;kwargs...) -> typeof(nv)

Return a copy of a concrete `NamedVector` with some of the field values replaced by the keyword arguments.
"""
Base.replace(nv::NamedVector;kwargs...)=replace(efficientoperations,nv;kwargs...)

"""
    map(f,nvs::NV...) where NV<:NamedVector -> NV

Apply function `f` elementwise on the input named vectors.
"""
@generated function Base.map(f,nvs::NV...) where NV<:NamedVector
    exprs=Vector{Expr}(undef,NV|>fieldcount)
    for i=1:length(exprs)
        tmp=[:(nvs[$j][$i]) for j=1:length(nvs)]
        exprs[i]=:(f($(tmp...)))
    end
    return :(($NV)($(exprs...)))
end

"""
    @namedvector structdef::Expr

Decorate a "raw" struct to be a subtype of `NamedVector`. Here, "raw" means that the input struct has no explicit supertype and no inner constructors.
"""
macro namedvector(structdef::Expr)
    tf=TypeFactory(structdef)
    @assert tf.supertype==Inference(:Any) "@namedvector error: no explicit supertype except `Any` is allowed."
    @assert length(tf.constructors)==0 "@namedvector error: no inner constructor is allowed."
    tf.supertype=Inference(:NamedVector)
    paramnames=tuple((param.name for param in tf.params)...)
    fieldnames=tuple((field.name for field in tf.fields)...)
    fldnm=FunctionFactory(name=:(Base.fieldnames))
    addargs!(fldnm,:(::Type{<:$(tf.name)}))
    extendbody!(fldnm,Expr(:tuple,QuoteNode.(fieldnames)...))
    structdef=tf(MixEscaped(Escaped(tf.name),UnEscaped(paramnames...,:NamedVector)))
    fldnmdef=fldnm(MixEscaped(Escaped(:tuple)))
    return Expr(:block,:(Base.@__doc__($structdef)),fldnmdef)
end

"""
    HomoNamedVector{T}

Abstract type for all homogeneous named vectors.
"""
abstract type HomoNamedVector{T} <: NamedVector end

"""
    eltype(::Type{<:HomoNamedVector{T}}) where T
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
