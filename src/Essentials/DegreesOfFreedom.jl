module DegreesOfFreedom

using Printf: @printf,@sprintf
using StaticArrays: SVector
using LaTeXStrings: latexstring
using ..Spatials: PID,AbstractBond
using ...Interfaces: rank,dimension
using ...Prerequisites: Float,decimaltostr,rawtype
using ...Prerequisites.TypeTraits: efficientoperations
using ...Prerequisites.CompositeStructures: CompositeDict
using ...Mathematics.VectorSpaces: VectorSpace
using ...Mathematics.AlgebraOverFields: SimpleID,ID,Element,Elements

import ..Spatials: pidtype,rcoord,icoord
import ...Interfaces: reset!,update!,sequence
import ...Mathematics.AlgebraOverFields: idpropertynames

export IID,Index,pid,iidtype,iid
export IndexToTuple,DirectIndexToTuple,directindextotuple,FilteredAttributes
export Internal,IDFConfig,Table
export LaTeX,OID,Operator,Operators,isHermitian,twist
export coordpresent,coordabsent
export latexformat,script,oidtype,otype
export Boundary

"""
    IID

The id of an internal degree of freedom.
"""
abstract type IID <: SimpleID end

"""
    Internal

The whole internal degrees of freedom at a single point.
"""
abstract type Internal{I<:IID} <: VectorSpace{I} end

"""
    show(io::IO,i::Internal)

Show an internal.
"""
Base.show(io::IO,i::Internal)=@printf io "%s(%s)" i|>typeof|>nameof join(("$name=$(getfield(i,name))" for name in i|>typeof|>fieldnames),",")

"""
    Index{P,I}

The complete index of a degree of freedom, which consist of the spatial part and the internal part.
"""
abstract type Index{P<:PID,I<:IID} <: SimpleID end

"""
    (INDEX::Type{<:Index})(pid::PID,iid::IID) -> INDEX

Get the corresponding index from a pid and an iid.
"""
@generated function (INDEX::Type{<:Index})(pid::PID,iid::IID)
    INDEX<:UnionAll ? :(INDEX(convert(Tuple,pid)...,convert(Tuple,iid)...)) : :(rawtype(INDEX)(convert(Tuple,pid)...,convert(Tuple,iid)...))
end

"""
    pidtype(index::Index)
    pidtype(::Type{<:Index{P}}) where {P<:PID}

Get the type of the spatial part of an index.
"""
pidtype(index::Index)=index|>typeof|>pidtype
pidtype(::Type{<:Index{P}}) where {P<:PID}=P

"""
    pid(index::Index) -> PID

Get the spatial part of an index.
"""
@generated function pid(index::Index)
    exprs=[:(getfield(index,$i)) for i=1:fieldcount(index|>pidtype)]
    return :(rawtype(pidtype(index))($(exprs...)))
end

"""
    iidtype(index::Index)
    iidtype(::Type{<:Index{<:PID,I}}) where {I<:IID}

Get the type of the internal part of an index.
"""
iidtype(index::Index)=index|>typeof|>iidtype
iidtype(::Type{<:Index{<:PID,I}}) where {I<:IID}=I

"""
    iid(index::Index) -> IID

Get the internal part of an index.
"""
@generated function iid(index::Index)
    exprs=[:(getfield(index,$i)) for i=fieldcount(index|>pidtype)+1:fieldcount(index)]
    return :(rawtype(iidtype(index))($(exprs...)))
end

"""
    adjoint(index::Index) -> typeof(index)

Get the adjoint of an index.
"""
Base.adjoint(index::Index)=typeof(index)(index|>pid,index|>iid|>adjoint)

"""
    union(::Type{P},::Type{I}) where {P<:PID,I<:IID}

Combine a concrete `PID` type and a concrete `IID` type to a concrete `Index` type.
"""
Base.union(::Type{P},::Type{I}) where {P<:PID,I<:IID}=Index{P,I}

"""
    IndexToTuple

The rules for converting an index to a tuple.

As a function, every instance should accept only one positional argument, i.e. the index to be converted to a tuple.
"""
abstract type IndexToTuple <:Function end
Base.:(==)(itt1::T,itt2::T) where T<:IndexToTuple = ==(efficientoperations,itt1,itt2)
Base.isequal(itt1::T,itt2::T) where T<:IndexToTuple = isequal(efficientoperations,itt1,itt2)

"""
    DirectIndexToTuple

Direct index to tuple.
"""
struct DirectIndexToTuple <: IndexToTuple end

"""
    (indextotuple::DirectIndexToTuple)(index::Index) -> Tuple

Convert an index to tuple directly.
"""
(indextotuple::DirectIndexToTuple)(index::Index)=convert(Tuple,index)

"""
    directindextotuple

Indicate that the conversion from an index to a tuple is direct.
"""
const directindextotuple=DirectIndexToTuple()

"""
    FilteredAttributes(::Type{I}) where I<:Index

A method that converts an arbitary index to a tuple, by iterating over the selected attributes in a specific order.
"""
struct FilteredAttributes{N} <: IndexToTuple
    attributes::NTuple{N,Symbol}
end
FilteredAttributes(attrs::Symbol...)=FilteredAttributes(attrs)
FilteredAttributes(::Type{I}) where I<:Index=FilteredAttributes(I|>fieldnames)

"""
    length(indextotuple::FilteredAttributes) -> Int
    length(::Type{<:FilteredAttributes{N}}) where N -> Int

Get the length of the filtered attributes.
"""
Base.length(indextotuple::FilteredAttributes)=indextotuple|>typeof|>length
Base.length(::Type{<:FilteredAttributes{N}}) where N=N

"""
    (indextotuple::FilteredAttributes)(index::Index) -> Tuple

Convert an index to tuple by the "filtered attributes" method.
"""
@generated function (indextotuple::FilteredAttributes)(index::Index)
    exprs=[:(getfield(index,indextotuple.attributes[$i])) for i=1:length(indextotuple)]
    return Expr(:tuple,exprs...)
end

"""
    filter(f::Function,indextotuple::FilteredAttributes)

Filter the attributes of a "filtered attributes" method.
"""
Base.filter(f::Function,indextotuple::FilteredAttributes)=FilteredAttributes(Tuple(attr for attr in indextotuple.attributes if f(attr)))

"""
    IDFConfig{I}(map::Function,pids::Union{AbstractVector{<:PID},Tuple{}}=()) where I<:Internal

Configuration of the internal degrees of freedom at a lattice.

`map` maps a `PID` to an `Internal`.
"""
struct IDFConfig{I<:Internal,M<:Function,P<:PID} <: CompositeDict{P,I}
    map::M
    contents::Dict{P,I}
end
function IDFConfig{I}(map::Function,pids::Union{AbstractVector{<:PID},Tuple{}}=()) where I<:Internal
    contents=Dict{pids|>eltype,I}()
    for pid in pids
        contents[pid]=map(pid)
    end
    IDFConfig(map,contents)
end

"""
    reset!(config::IDFConfig,pids) -> IDFConfig

Reset the idfconfig with new pids.
"""
function reset!(config::IDFConfig,pids)
    empty!(config)
    for pid in pids
        config[pid]=config.map(pid)
    end
    config
end

"""
    Table{I}(by::IndexToTuple) where I<:Index

Index-sequence table.
"""
struct Table{I<:Index,B<:IndexToTuple} <: CompositeDict{I,Int}
    by::B
    contents::Dict{I,Int}
end
Table{I}(by::IndexToTuple) where I<:Index=Table(by,Dict{I,Int}())

"""
    Table(indices::AbstractVector{<:Index},by::IndexToTuple=directindextotuple) -> Table

Convert an sequence of indices to the corresponding index-sequence table.

The input indices will be converted to tuples by the `by` function with the duplicates removed. The resulting unique tuples are sorted, which determines the sequence of the input indices. Note that two indices have the same sequence if their converted tupels are equal to each other.
"""
function Table(indices::AbstractVector{<:Index},by::IndexToTuple=directindextotuple)
    tuples=[by(index) for index in indices]
    permutation=sortperm(tuples,alg=Base.Sort.QuickSort)
    result=Table{indices|>eltype}(by)
    count=1
    for i=1:length(tuples)
        i>1 && tuples[permutation[i]]!=tuples[permutation[i-1]] && (count+=1)
        result[indices[permutation[i]]]=count
    end
    result
end

"""
    Table(config::IDFConfig,by::IndexToTuple=directindextotuple) -> Table

Get the index-sequence table of the whole internal Hilbert spaces at a lattice.
"""
function Table(config::IDFConfig,by::IndexToTuple=directindextotuple)
    result=union(config|>keytype,config|>valtype|>eltype)[]
    for (pid,internal) in config
        for iid in internal
            push!(result,(result|>eltype)(pid,iid))
        end
    end
    Table(result,by)
end

"""
    union(tables::Table...) -> Table

Unite several index-sequence tables.

See [`Table`](@ref) for more details.
"""
function Base.union(tables::Table...)
    @assert mapreduce(table->table.by,==,tables)
    indices=(tables|>eltype|>keytype)[]
    for table in tables
        for index in keys(table)
            push!(indices,index)
        end
    end
    Table(indices,tables[1].by)
end

"""
    reverse(table::Table) -> Dict{Int,Set{<:Index}}

Convert an index-sequence table to a sequence-indices table.

Since different indices may correspond to the same sequence, the reverse is a one-to-many map.
"""
function Base.reverse(table::Table)
    result=Dict{Int,Set{table|>keytype}}()
    for (index,seq) in table
        haskey(result,seq) || (result[seq]=Set{table|>keytype}())
        push!(result[seq],index)
    end
    result
end

"""
    reset!(table::Table,indices::AbstractVector{<:Index}) -> Table
    reset!(table::Table,config::IDFConfig) -> Table

Reset a table.
"""
function reset!(table::Table,indices::AbstractVector{<:Index})
    empty!(table)
    tuples=[table.by(index) for index in indices]
    permutation=sortperm(tuples,alg=Base.Sort.QuickSort)
    count=1
    for i=1:length(tuples)
        i>1 && tuples[permutation[i]]!=tuples[permutation[i-1]] && (count+=1)
        table[indices[permutation[i]]]=count
    end
    table
end
function reset!(table::Table,config::IDFConfig)
    indices=union(config|>keytype,config|>valtype|>eltype)[]
    for (pid,internal) in config
        for iid in internal
            push!(indices,(indices|>eltype)(pid,iid))
        end
    end
    reset!(table,indices)
end

"""
    LaTeX{SP,SB}(body) where {SP,SB}

LaTeX string representation.
"""
struct LaTeX{SP,SB,B}
    body::B
    function LaTeX{SP,SB}(body=nothing) where {SP,SB}
        @assert isa(SP,Tuple{Vararg{Symbol}}) && isa(SB,Tuple{Vararg{Symbol}}) "LaTeX error: SP and SB must be tuple of symbols."
        new{SP,SB,typeof(body)}(body)
    end
end
latexsuperscript(::Type{<:LaTeX{SP}}) where SP=SP
latexsubscript(::Type{<:LaTeX{SP,SB} where SP}) where SB=SB

@generated function oidcoord(vector::SVector{N,Float}) where N
    exprs=[:(vector[$i]===-0.0 ? 0.0 : vector[$i]) for i=1:N]
    return :(SVector($(exprs...)))
end

"""
    OID(index::Index,::Nothing,::Nothing,seq::Union{Nothing,Int})
    OID(index::Index,rcoord::SVector{N,Float},icoord::SVector{N,Float},seq::Union{Nothing,Int}) where N
    OID(index::Index,rcoord::Vector{Float},icoord::Vector{Float},seq::Union{Nothing,Int})
    OID(index::Index;rcoord::Union{Nothing,SVector,Vector{Float}}=nothing,icoord::Union{Nothing,SVector,Vector{Float}}=nothing,seq::Union{Nothing,Int}=nothing)

Operator id.
"""
struct OID{I<:Index,RC<:Union{Nothing,SVector},IC<:Union{Nothing,SVector},S<:Union{Nothing,Int}} <: SimpleID
    index::I
    rcoord::RC
    icoord::IC
    seq::S
    OID(index::Index,::Nothing,::Nothing,seq::Union{Nothing,Int})=new{typeof(index),Nothing,Nothing,typeof(seq)}(index,nothing,nothing,seq)
    function OID(index::Index,rcoord::SVector{N,Float},icoord::SVector{N,Float},seq::Union{Nothing,Int}) where N
        new{typeof(index),SVector{N,Float},SVector{N,Float},typeof(seq)}(index,oidcoord(rcoord),oidcoord(icoord),seq)
    end
end
OID(index::Index,rcoord::Vector{Float},icoord::Vector{Float},seq::Union{Nothing,Int})=OID(index,SVector{length(rcoord)}(rcoord),SVector{length(icoord)}(icoord),seq)
function OID(index::Index;rcoord::Union{Nothing,SVector,Vector{Float}}=nothing,icoord::Union{Nothing,SVector,Vector{Float}}=nothing,seq::Union{Nothing,Int}=nothing)
    OID(index,rcoord,icoord,seq)
end
totuple(::Nothing)=nothing
@generated totuple(v::SVector{N}) where N=Expr(:tuple,[:(v[$i]) for i=1:N]...)
Base.hash(oid::OID{<:Index},h::UInt)=hash((oid.index,totuple(oid.rcoord)),h)
Base.fieldnames(::Type{<:OID})=(:index,:rcoord,:icoord,:seq)
idpropertynames(::Type{<:ID{OID}})=(:indexes,:rcoords,:icoords,:seqs)

"""
    show(io::IO,oid::OID)

Show an operator id.
"""
function Base.show(io::IO,oid::OID)
    @printf io "OID(%s" oid.index
    oid.rcoord===nothing ? (@printf io ",:") : (@printf io ",[%s]" join(oid.rcoord,","))
    oid.icoord===nothing ? (@printf io ",:") : (@printf io ",[%s]" join(oid.icoord,","))
    @printf io ",%s)" oid.seq===nothing ? ":" : oid.seq
end

"""
    repr(oid::OID,l::LaTeX) -> String

LaTeX string representation of an oid.
"""
Base.repr(oid::OID,l::LaTeX)=@sprintf "%s^{%s}_{%s}" script(oid,l,Val(:B)) join(script(oid,l,Val(:SP)),',') join(script(oid,l,Val(:SB)),',')

"""
    adjoint(oid::OID) -> typeof(oid)
    adjoint(oid::ID{OID,N}) where N -> typeof(oid)

Get the adjoint of an operator id.
"""
Base.adjoint(oid::OID)=OID(oid.index',oid.rcoord,oid.icoord,oid.seq)
@generated Base.adjoint(oid::ID{OID,N}) where N=Expr(:call,:ID,[:(oid[$i]') for i=N:-1:1]...)

"""
    isHermitian(oid::ID{OID,N}) where N -> Bool

Judge whether an operator id is Hermitian.
"""
function isHermitian(oid::ID{OID,N}) where N
    for i=1:((N+1)÷2)
        oid[i]'==oid[N+1-i] || return false
    end
    return true
end

"""
    script(oid::OID,::Val{attr}) where attr -> String

Get the `:rcoord/:icoord` script of an oid.
"""
@generated function script(oid::OID,::Val{attr}) where attr
    @assert attr in (:rcoord,:icoord) "script error: not supported attr($attr)."
    fieldtype(oid,attr)<:Nothing ? 'N' : (attr=QuoteNode(attr);:("[$(join(getfield(oid,$attr),','))]"))
end

"""
    script(oid::OID,l::LaTeX,::Val{:B}) -> Any
    script(oid::OID,l::LaTeX,::Val{:SP}) -> Tuple
    script(oid::OID,l::LaTeX,::Val{:SB}) -> Tuple

Get the body/superscript/subscript of the latex string representation of an oid.
"""
@generated script(oid::OID,l::LaTeX,::Val{:B})=fieldtype(l,:body)<:Nothing ? :(oid.index.scope) : :(l.body)
@generated script(oid::OID,l::LaTeX,::Val{:SP})=Expr(:tuple,[:(script(oid,Val($super))) for super in QuoteNode.(l|>latexsuperscript)]...)
@generated script(oid::OID,l::LaTeX,::Val{:SB})=Expr(:tuple,[:(script(oid,Val($sub))) for sub in QuoteNode.(l|>latexsubscript)]...)

"""
    coordpresent

Indicate that the `:icoord` and `:rcoord` attributes in an oid should not be nothing.
"""
const coordpresent=Val(true)

"""
    coordabsent

Indicate that the `:icoord` and `:rcoord` attributes in an oid should be nothing
"""
const coordabsent=Val(false)

"""
    oidtype(I::Type{<:IID},B::Type{<:AbstractBond},::Type{Nothing},::Val{true})
    oidtype(I::Type{<:IID},B::Type{<:AbstractBond},::Type{<:Table},::Val{true})
    oidtype(I::Type{<:IID},B::Type{<:AbstractBond},::Type{Nothing},::Val{false})
    oidtype(I::Type{<:IID},B::Type{<:AbstractBond},::Type{<:Table},::Val{false})

Get the compatible oid type.
"""
oidtype(I::Type{<:IID},B::Type{<:AbstractBond},::Type{Nothing},::Val{true})=OID{union(B|>pidtype,I),SVector{B|>dimension,Float},SVector{B|>dimension,Float},Nothing}
oidtype(I::Type{<:IID},B::Type{<:AbstractBond},::Type{<:Table},::Val{true})=OID{union(B|>pidtype,I),SVector{B|>dimension,Float},SVector{B|>dimension,Float},Int}
oidtype(I::Type{<:IID},B::Type{<:AbstractBond},::Type{Nothing},::Val{false})=OID{union(B|>pidtype,I),Nothing,Nothing,Nothing}
oidtype(I::Type{<:IID},B::Type{<:AbstractBond},::Type{<:Table},::Val{false})=OID{union(B|>pidtype,I),Nothing,Nothing,Int}

"""
    angle(id::ID{OID},vectors::AbstractVector{<:AbstractVector{Float}},values::AbstractVector{Float}) -> Float

Get the total twist phase of an id.
"""
@generated function Base.angle(id::ID{OID},vectors::AbstractVector{<:AbstractVector{Float}},values::AbstractVector{Float})
    Expr(:call,:+,[:(angle(id[$i],vectors,values)) for i=1:rank(id)]...)
end

"""
    Operator{V<:Number,I<:ID{OID}} <: Element{V,I}

Abstract type for an operator.
"""
abstract type Operator{V<:Number,I<:ID{OID}} <: Element{V,I} end
function (O::Type{<:Operator})( value,
                                indexes::NTuple{N,Index};
                                rcoords::Union{Nothing,NTuple{N,SVector{M,Float}}}=nothing,
                                icoords::Union{Nothing,NTuple{N,SVector{M,Float}}}=nothing,
                                seqs::Union{Nothing,NTuple{N,Int}}=nothing
                                ) where {N,M}
    rcoords===nothing && (rcoords=ntuple(i->nothing,N))
    icoords===nothing && (icoords=ntuple(i->nothing,N))
    seqs===nothing && (seqs=ntuple(i->nothing,N))
    return O(value,ID(OID,indexes,rcoords,icoords,seqs))
end

const latexformats=Dict{Symbol,LaTeX}()
"""
    latexformat(T::Type{<:Operator}) -> LaTeX
    latexformat(T::Type{<:Operator},l::LaTeX) -> LaTeX

Get/Set the latex format for a subtype of `Operator`.
"""
latexformat(T::Type{<:Operator})=latexformats[nameof(T)]
latexformat(T::Type{<:Operator},l::LaTeX)=latexformats[nameof(T)]=l

"""
    show(io::IO,opt::Operator)
    show(io::IO,::MIME"text/latex",opt::Operator)

Show an operator.
"""
Base.show(io::IO,opt::Operator)=@printf io "%s(%s,%s)" nameof(typeof(opt)) decimaltostr(opt.value) opt.id
Base.show(io::IO,::MIME"text/latex",opt::Operator)=show(io,MIME"text/latex"(),latexstring(repr(opt)))

"""
    adjoint(opt::Operator{N}) where N -> Operator

Get the adjoint of an operator.
"""
Base.adjoint(opt::Operator)=rawtype(typeof(opt))(opt.value',opt.id')

"""
    repr(opt::Operator,l::Union{LaTeX,Nothing}=nothing) -> String

Get the latex string representation of an operator.
"""
Base.repr(opt::Operator,::Nothing=nothing)=repr(opt,latexformat(typeof(opt)))
function Base.repr(opt::Operator,l::LaTeX)
    rank(opt)==0 && return replace(valuetolatextext(opt.value)," "=>"")
    return @sprintf "%s%s" valuetostr(opt.value) join(NTuple{rank(opt),String}(repr(opt.id[i],l) for i=1:rank(opt)),"")
end
function valuetostr(v)
    v==1 && return ""
    v==-1 && return "-"
    result=valuetolatextext(v)
    if occursin(" ",result) || (isa(v,Complex) && real(v)≠0 && imag(v)≠0)
        bra,ket=occursin("(",result) ? ("[","]") : ("(",")")
        result=@sprintf "%s%s%s" bra replace(result," "=>"") ket
    end
    return result
end
valuetolatextext(value::Union{Real,Complex})=decimaltostr(value)
function valuetolatextext(value)
    io=IOBuffer()
    showable(MIME"text/latex"(),value) ? show(IOContext(io,:limit=>false),MIME"text/latex"(),value) : show(IOContext(io,:limit=>false),MIME"text/plain"(),value)
    return replace(replace(String(take!(io)),"\\begin{equation*}"=>""),"\\end{equation*}"=>"")
end

"""
    isHermitian(opt::Operator) -> Bool

Judge whether an operator is Hermitian.
"""
isHermitian(opt::Operator)=isa(opt.value,Real) && isHermitian(opt.id)

"""
    rcoord(opt::Operator) -> SVector

Get the whole rcoord of an operator.
"""
@generated function rcoord(opt::Operator)
    rank(opt)==1 && return :(opt.id[1].rcoord)
    rank(opt)==2 && return :(opt.id[1].rcoord-opt.id[2].rcoord)
    error("rcoord error: not supported rank($(rank(opt))) of $(nameof(opt)).")
end

"""
    icoord(opt::Operator) -> SVector

Get the whole icoord of an operator.
"""
@generated function icoord(opt::Operator)
    rank(opt)==1 && return :(opt.id[1].icoord)
    rank(opt)==2 && return :(opt.id[1].icoord-opt.id[2].icoord)
    error("icoord error: not supported rank($(rank(opt))) of $(nameof(opt)).")
end

"""
    sequence(opt::Operator,table=nothing) -> NTuple{rank(opt),Int}

Get the sequence of the oids of an operator according to a table.
"""
@generated function sequence(opt::Operator,table=nothing)
    table<:Nothing && return :(opt.id.seqs::NTuple{rank(opt),Int})
    return Expr(:tuple,[:(get(table,opt.id[$i].index,nothing)) for i=1:rank(opt)]...)
end

"""
    otype

Get the compatible operator type from a term type, a bond type and a table type.
"""
function otype end

"""
    Operators(opts::Operator...)

A set of operators.

Type alias for `Elements{<:ID{OID},<:Operator}`.
"""
const Operators{I<:ID{OID},O<:Operator}=Elements{I,O}
Operators(opts::Operator...)=Elements(opts...)

"""
    adjoint(opts::Operators) -> Operators

Get the adjoint of a set of operators.
"""
function Base.adjoint(opts::Operators)
    result=Operators{opts|>keytype,opts|>valtype}()
    for opt in values(opts)
        nopt=opt|>adjoint
        result[nopt.id]=nopt
    end
    return result
end

"""
    repr(opts::Operators,l::Union{LaTeX,Nothing}=nothing) -> String

Get the latex string representation of a set of operators.
"""
function Base.repr(opts::Operators,l::Union{LaTeX,Nothing}=nothing)
    result=String[]
    for (i,opt) in enumerate(values(opts))
        rep=repr(opt,l)
        i>1 && rep[1]!='-' && push!(result,"+")
        push!(result,rep)
    end
    return join(result,"")
end

"""
    show(io::IO,::MIME"text/latex",opts::Operators)

Show latex formed operators.
"""
Base.show(io::IO,::MIME"text/latex",opts::Operators)=show(io,MIME"text/latex"(),latexstring(repr(opts)))

"""
    summary(io::IO,opts::Operators)

Print a brief description of a set of operators to an io.
"""
Base.summary(io::IO,opts::Operators)=@printf io "Operators{%s}" valtype(opts)

"""
    isHermitian(opts::Operators) -> Bool

Judge whether a set of operators as a whole is Hermitian.
"""
isHermitian(opts::Operators)=opts==opts'

"""
    twist(operator::Operator,vectors::AbstractVector{<:AbstractVector{Float}},values::AbstractVector{Float}) -> Operator

Twist an operator.
"""
function twist(operator::Operator,vectors::AbstractVector{<:AbstractVector{Float}},values::AbstractVector{Float})
    replace(operator,value=operator.value*exp(1im*angle(operator.id,vectors,values)))
end

"""
    Boundary{Names}(values::AbstractVector{Float},vectors::AbstractVector{<:AbstractVector{Float}}) where Names
    Boundary()

Boundary twist of operators.
"""
struct Boundary{Names,V<:AbstractVector{Float}}
    values::Vector{Float}
    vectors::Vector{V}
    function Boundary{Names}(values::AbstractVector{Float},vectors::AbstractVector{<:AbstractVector{Float}}) where Names
        @assert length(Names)==length(values)==length(vectors) "Boundary error: dismatched names, values and vectors."
        new{Names,eltype(vectors)}(convert(Vector{Float},values),vectors)
    end
end
const boundaryemptyvalues=Float[]
const boundaryemptyvectors=SVector{0,Float}[]
Boundary()=Boundary{()}(boundaryemptyvalues,boundaryemptyvectors)

"""
    ==(bound1::Boundary,bound2::Boundary) -> Bool
    isequal(bound1::Boundary,bound2::Boundary) -> Bool
"""
Base.:(==)(bound1::Boundary,bound2::Boundary) = ==(efficientoperations,bound1,bound2)
Base.isequal(bound1::Boundary,bound2::Boundary)=isequal(efficientoperations,bound1,bound2)

"""
    (bound::Boundary)(operator::Operator) -> Operator
    (bound::Boundary{()})(operator::Operator) -> Operator

Get the boundary twisted operator.
"""
(bound::Boundary)(operator::Operator)=twist(operator,bound.vectors,bound.values)
(bound::Boundary{()})(operator::Operator)=operator

"""
    angle(bound::Boundary,operator::Operator) -> Float
    angle(bound::Boundary{()},operator::Operator) -> Int

Get the boundary twist phase of an operator.
"""
Base.angle(bound::Boundary,operator::Operator)=angle(operator.id,bound.vectors,bound.values)
Base.angle(bound::Boundary{()},operator::Operator)=0

"""
    update!(bound::Boundary{Names},args...;kwargs...) where Names -> Boundary

Update the values of the boundary twisted phase.
"""
@generated function update!(bound::Boundary{Names},args...;kwargs...) where Names
    return Expr(:block,[:(bound.values[$i]=get(kwargs,$name,bound.values[$i])) for (i,name) in enumerate(QuoteNode.(Names))]...,:(bound))
end

end #module
