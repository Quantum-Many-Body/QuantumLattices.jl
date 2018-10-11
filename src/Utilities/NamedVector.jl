module NamedVector

import Printf: @printf

export AbstractNamedVector
export @namedvector

"""
    AbstractNamedVector{T,A<:AbstractVector{T}}

Abstract type for all concrete named vectors. To subtype it, please note:
1. The concrete types must have the field `values::A`, which is used to store the values by design. A recommended template for the subtype is
    ```
    struct YourNamedVector{T,A} <: AbstractNamedVector{T,A}
        values::A
    end
    ```
2. The concrete types must implement their own `Base.fieldnames`, which defines the names of the type. A recommended template for this method is
    ```
    Base.fieldnames(::Type{YourNamedVector},private=false)=private ? tuple(YourNames...,:values) : YourNames
    ```
3. Arithmetic operations, such as `+`, `-`, `*`, `\\`, `%`, `÷`, etc. are **NOT** supported. Thus users have to overload these operations themselves if needed.
"""
abstract type AbstractNamedVector{T,A<:AbstractVector{T}} end

"""
    getproperty(nv::AbstractNamedVector,key::Symbol)

Get the value by the `.` syntax.
"""
function Base.getproperty(nv::AbstractNamedVector,key::Symbol)
    if key==:values
        result=getfield(nv,key)
    elseif key ∈ typeof(nv)|>fieldnames
        result=getfield(nv,:values)[findfirst(isequal(key),typeof(nv)|>fieldnames)]
    else
        error("$(string(typeof(nv))) getproperty error: :$key not available.")
    end
    result
end

"""
    setproperty!(nv::AbstractNamedVector,key::Symbol,value)

Set the value by the `.` syntax. Note the value of the field `:values` as a whole can not be changed directly by design. Instead, the elements of `:values` can be changed.
"""
function Base.setproperty!(nv::AbstractNamedVector,key::Symbol,value)
    @assert key ∈ typeof(nv)|>fieldnames "setproperty! error: set the value of :$key not available."
    getfield(nv,:values)[findfirst(isequal(key),typeof(nv)|>fieldnames)]=convert(nv|>typeof|>eltype,value)
end

"""
    getindex(nv::AbstractNamedVector,index::Int)

Get the value by the `[]` syntax.
"""
Base.getindex(nv::AbstractNamedVector,index::Int)=getfield(nv,:values)[index]

"""
    setindex!(nv::AbstractNamedVector,value,index::Int)

Set the value by the `[]` syntax.
"""
Base.setindex!(nv::AbstractNamedVector,value,index::Int)=(getfield(nv,:values)[index]=convert(nv|>typeof|>eltype,value))

"""
    ==(nv1::AbstractNamedVector,nv2::AbstractNamedVector)

Overloaded `==` operator.
"""
Base.:(==)(nv1::AbstractNamedVector,nv2::AbstractNamedVector)=nv1|>typeof|>fieldnames==nv2|>typeof|>fieldnames && getfield(nv1,:values)==getfield(nv2,:values)

"""
    show(io::IO,nv::AbstractNamedVector)

Show a concrete `AbstractNamedVector`.
"""
Base.show(io::IO,nv::AbstractNamedVector)=@printf io "%s(%s)" Base.typename(typeof(nv)) join(getfield(nv,:values),',')

"""
    hash(nv::AbstractNamedVector,h::UInt)

Hash a concrete `AbstractNamedVector`.
"""
Base.hash(nv::AbstractNamedVector,h::UInt)=hash(getfield(nv,:values),h)

"""
    length(nv::AbstractNamedVector)

Get the length of a concrete `AbstractNamedVector`.
"""
Base.length(nv::AbstractNamedVector)=length(getfield(nv,:values))

"""
    eltype(::Type{NV},choice::Int=1) where NV<:AbstractNamedVector{T,A} where {T,A}

Get the type parameter of a concrete `AbstractNamedVector`.
"""
Base.eltype(::Type{NV},choice::Int=1) where NV<:AbstractNamedVector{T,A} where {T,A}=(0<choice<3 || error("eltype error: choice($choice) must be 1 or 2.");choice==1 ? T : A)

"""
    zero(nv::AbstractNamedVector)

Get a concrete `AbstractNamedVector` with all values being zero.
"""
function Base.zero(nv::AbstractNamedVector)
    result=deepcopy(nv)
    for (i,value) in enumerate(getfield(nv,:values))
        result[i]=zero(value)
    end
    result
end

"""
    iterate(nv::AbstractNamedVector)
    iterate(nv::AbstractNamedVector,state)

Iterate over the values of a concrete `AbstractNamedVector`.
"""
Base.iterate

Base.iterate(nv::AbstractNamedVector)=iterate(getfield(nv,:values))
Base.iterate(nv::AbstractNamedVector,state)=iterate(getfield(nv,:values),state)

"""
    keys(nv::AbstractNamedVector)

Iterate over the names.
"""
Base.keys(nv::AbstractNamedVector)=(name for name in nv|>typeof|>fieldnames)

"""
    values(nv::AbstractNamedVector)

Iterate over the values.
"""
Base.values(nv::AbstractNamedVector)=(value for value in getfield(nv,:values))

"""
    pairs(nv::AbstractNamedVector)

Iterate over the `name=>value` pairs.
"""
Base.pairs(nv::AbstractNamedVector)=Base.Generator(=>,keys(nv),values(nv))

"""
    replace(nv::AbstractNamedVector;kwargs...)

Return a copy of a concrete `AbstractNamedVector` with some of the filed values replaced by the keyword arguments.
"""
function Base.replace(nv::AbstractNamedVector;kwargs...)
    result=deepcopy(nv)
    for (i,(name,value)) in enumerate(zip(nv|>typeof|>fieldnames,getfield(nv,:values)))
        result[i]=get(kwargs,name,value)
    end
    result
end

"""
    @namedvector typename fieldnames dtype::Union{Expr,Symbol}=:nothing atype::Union{Expr,Symbol}=:nothing supertypename=:AbstractNamedVector

Construct a concrete named vector with the type name being `typename`, fieldnames specified by `fieldnames`, and optionally, the type parameters specified by `dtype` and `atype`, and the supertype specified by `supertypename`.
"""
macro namedvector(typename,fieldnames,dtype::Union{Expr,Symbol}=:nothing,atype::Union{Expr,Symbol}=:nothing,supertypename=:AbstractNamedVector)
    typename,supertypename=Symbol(typename),Symbol(supertypename)
    dsp=supertypename==:AbstractNamedVector
    fieldnames=tuple(eval(fieldnames)...)
    @assert all(isa(name,Symbol) && name!=:values for name in fieldnames) "namedvector error: every field name should be a `Symbol` but not be `:values`."
    isa(dtype,Expr) && (@assert (dtype.head==:(<:) && dtype.args|>length==1) "namedvector error: not supported `dtype`.")
    isa(atype,Expr) && (@assert (atype.head==:(<:) && atype.args|>length==1) "namedvector error: not supported `atype`.")
    dname,dscope=isa(dtype,Expr) ? (:T,dtype.args[1]) : (dtype==:nothing ? (:T,:Any) : (dtype,:concrete))
    aname,ascope=isa(atype,Expr) ? (:A,atype.args[1]) : (atype==:nothing ? (:A,:AbstractVector) : (atype,:concrete))
    if dscope==:concrete && ascope==:concrete
        return quote
            struct $(esc(typename)) <: ($dsp ? AbstractNamedVector : $(esc(supertypename))){$(esc(dname)),$(esc(aname)){$(esc(dname))}}
                values::$(esc(aname)){$(esc(dname))}
                $(esc(typename))(values::$(esc(aname)){$(esc(dname))})=new(values)
                $(esc(typename))(values::$(esc(dname))...)=$(esc(typename))(convert($(esc(aname)),collect(values)))
            end
            Base.fieldnames(::Type{$(esc(typename))},private=false)=private ? tuple(($fieldnames)...,:values) : $fieldnames
        end
    elseif dscope==:concrete
        return quote
            struct $(esc(typename)){$aname<:$(esc(ascope)){$(esc(dname))}} <: ($dsp ? AbstractNamedVector : $(esc(supertypename))){$(esc(dname)),$aname}
                values::$aname
                $(esc(typename))(values::$(esc(ascope)){$(esc(dname))})=new{values|>typeof}(values)
                $(esc(typename))(values::$(esc(dname))...)=$(esc(typename))(collect(values))
            end
            $(esc(typename)){$aname}(values::$(esc(ascope)){$(esc(dname))}) where $aname=$(esc(typename))(values)
            Base.fieldnames(::Type{<:$(esc(typename))},private=false)=private ? tuple(($fieldnames)...,:values) : $fieldnames
        end
    elseif ascope==:concrete
        return quote
            struct $(esc(typename)){$dname<:$(esc(dscope))} <: ($dsp ? AbstractNamedVector : $(esc(supertypename))){$dname,$(esc(aname)){$dname}}
                values::$(esc(aname)){$dname}
                $(esc(typename))(values::$(esc(aname)){<:$(esc(dscope))})=new{values|>eltype}(values)
                $(esc(typename))(values::$(esc(dscope))...)=$(esc(typename))(convert($(esc(aname)),collect(values)))
            end
            $(esc(typename)){$dname}(values::$(esc(aname)){<:$(esc(dscope))}) where $dname=$(esc(typename))(values)
            Base.fieldnames(::Type{<:$(esc(typename))},private=false)=private ? tuple(($fieldnames)...,:values) : $fieldnames
        end
    else
        return quote
            struct $(esc(typename)){$dname<:$(esc(dscope)),$aname<:$(esc(ascope)){$dname}} <: ($dsp ? AbstractNamedVector : $(esc(supertypename))){$dname,$aname}
                values::$aname
                $(esc(typename))(values::$(esc(ascope)){<:$(esc(dscope))})=new{values|>eltype,values|>typeof}(values)
                $(esc(typename))(values::$(esc(dscope))...)=$(esc(typename))(collect(values))
            end
            $(esc(typename)){$dname,$aname}(values::$(esc(ascope)){<:$(esc(dscope))}) where {$dname,$aname}=$(esc(typename))(values)
            Base.fieldnames(::Type{<:$(esc(typename))},private=false)=private ? tuple(($fieldnames)...,:values) : $fieldnames
        end
    end
end

end #module
