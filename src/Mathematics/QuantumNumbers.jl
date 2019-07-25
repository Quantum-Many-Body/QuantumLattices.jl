module QuantumNumbers

using Base.Iterators: Reverse,flatten,product,reverse
using Printf: @printf,@sprintf
using DataStructures: OrderedDict
using Random: MersenneTwister,seed!,shuffle!
using LinearAlgebra: norm
using ..Combinatorics: Combinations
using ...Prerequisites: Float
using ...Prerequisites.NamedVectors: HomoNamedVector

import ...Interfaces: ⊕,⊗,dimension,expand,permute,decompose,regularize!,regularize

export AbelianNumber,@abeliannumber,SQN,PQN,SPQN,Z2QN,periods
export qnscounts,qnsindptr,qnscompression,qnsexpansion,qnscontents,qnsindices,qnsbruteforce,qnsmontecarlo
export AbelianNumbers,SQNS,PQNS,SzPQNS,SPQNS,Z2QNS
export ukron,toordereddict

"Abstract type for all concrete quantum numbers for a single basis."
abstract type AbelianNumber <: HomoNamedVector{Float} end

"""
    dimension(::Type{<:AbelianNumber}) -> Int
    dimension(::AbelianNumber) -> Int

The dimension of the Hilbert space an `AbelianNumber` represents. Apparently, this is always 1.
"""
dimension(::Type{<:AbelianNumber})=1
dimension(::AbelianNumber)=1

"""
    regularize!(::Type{QN},array::AbstractVector{<:Real}) where QN<:AbelianNumber -> typeof(array)
    regularize!(::Type{QN},array::AbstractMatrix{<:Real}) where QN<:AbelianNumber -> typeof(array)

Regularize the elements of an array in place so that it can represent quantum numbers.
"""
function regularize!(::Type{QN},array::AbstractVector{<:Real}) where QN<:AbelianNumber
    @assert array|>length==QN|>length "regularize! error: not consistent shape of input array and $QN."
    for (i,period) in enumerate(QN|>periods)
        period==Inf || @inbounds begin
            remainder=array[i]%period
            array[i]=remainder<0 ? remainder+period : remainder
        end
    end
    return array
end
function regularize!(::Type{QN},array::AbstractMatrix{<:Real}) where QN<:AbelianNumber
    @assert size(array,1)==QN|>length "regularize! error: not consistent shape of input array and $QN."
    for (i,period) in enumerate(QN|>periods)
        period==Inf || @inbounds for j in 1:size(array,2)
            remainder=array[i,j]%period
            array[i,j]=remainder<0 ? remainder+period : remainder
        end
    end
    return array
end

"""
    regularize(::Type{QN},array::Union{AbstractVector{<:Real},AbstractMatrix{<:Real}}) where {QN<:AbelianNumber} -> typeof(array)

Regularize the elements of an array and return a copy that can represent quantum numbers.
"""
function regularize(::Type{QN},array::Union{AbstractVector{<:Real},AbstractMatrix{<:Real}}) where {QN<:AbelianNumber}
    result=copy(array)
    regularize!(QN,result)
    return result
end

"""
    @abeliannumber typename fieldnames fieldperiods

Construct a concrete `AbelianNumber` with the type name being `typename`, fieldnames specified by `fieldnames` and periods specified by `fieldperiods`.
"""
macro abeliannumber(typename,fieldnames,fieldperiods)
    typename=Symbol(typename)
    fieldnames=tuple(eval(fieldnames)...)
    fieldperiods=tuple(eval(fieldperiods)...)
    arguments=ntuple(i->Symbol(:v,i),length(fieldnames))
    @assert length(fieldnames)==length(fieldperiods) "abeliannumber error: number of fieldnames($(length(fieldnames))) and fieldperiods($(length(fieldperiods))) not equal."
    @assert all(isa(name,Symbol) for name in fieldnames) "abeliannumber error: all field names should be Symbol."
    @assert all(fieldperiods.>0) "abeliannumber error: all field periods should be positive."
    if all(fieldperiods.==Inf)
        title=Expr(:call,:($(esc(typename))),(:($arg::Float) for arg in arguments)...)
        body=Expr(:call,:new,arguments...)
    else
        title=Expr(:call,:($(esc(typename))),(:($arg::Float) for arg in arguments)...,Expr(:kw,:(regularize::Bool),:true))
        body=Expr(:call,:new,(p==Inf ? arg : :(regularize ? ($arg)%($p)+(($arg)%($p)<0 ? $p : 0) : $arg) for (p,arg) in zip(fieldperiods,arguments))...)
    end
    newtype=Expr(:struct,false,:($(esc(typename))<:AbelianNumber),Expr(:block,(:($field::Float) for field in fieldnames)...,Expr(:(=),title,body)))
    functions=Expr(:block,:(Base.fieldnames(::Type{<:$(esc(typename))})=$fieldnames),:(periods(::Type{<:$(esc(typename))})=$fieldperiods))
    return Expr(:block,:(global periods),:(Base.@__doc__($newtype)),functions)
end

"""
    SQN(Sz::Real)

The concrete `AbelianNumber` of a quantum system with spin z-component `Sz` conserved.
"""
@abeliannumber "SQN" (:Sz,) (Inf,)

"""
    PQN(N::Real)

The concrete `AbelianNumber` of a quantum system with particle number `N` conserved.
"""
@abeliannumber "PQN" (:N,) (Inf,)

"""
    SPQN(N::Real,Sz::Real)

The concrete `AbelianNumber` of a quantum system with both particle number `N` and spin z-component `Sz` conserved.
"""
@abeliannumber "SPQN" (:N,:Sz) (Inf,Inf)

"""
    Z2QN(N::Real)

The concrete `AbelianNumber` of a quantum system with a Z₂-like conserved quantity.
"""
@abeliannumber "Z2QN" (:N,) (2,)

abstract type QnsProtocol end
struct QnsIndptr <: QnsProtocol end
struct QnsCounts <: QnsProtocol end
"""
    qnsindptr

Indicate that methods with `AbelianNumbers` use the index pointer of the compressed contents.
"""
const qnsindptr=QnsIndptr()
"""
    qnscounts

Indicate that methods with `AbelianNumbers` use the count number of the compressed contents.
"""
const qnscounts=QnsCounts()

"""
    AbelianNumbers(form::Char,contents::Vector{<:AbelianNumber},counts::Vector{Int},::QnsCounts)
    AbelianNumbers(form::Char,contents::Vector{<:AbelianNumber},indptr::Vector{Int},::QnsIndptr)

The whole quantum numbers of the total bases of a Hilbert space.

The default constructors construct a `AbelianNumbers` from a vector of concrete quantum numbers and an vector containing their counts or indptr.
"""
struct AbelianNumbers{QN<:AbelianNumber}
    form::Char
    contents::Vector{QN}
    indptr::Vector{Int}
    function AbelianNumbers(form::Char,contents::Vector{<:AbelianNumber},counts::Vector{Int},::QnsCounts)
        @assert form|>uppercase ∈ ('G','U','C') "AbelianNumbers error: 'form'($form) is not 'G','U' or 'C'."
        @assert length(contents)==length(counts) "AbelianNumbers error: dismatch lengths of contents and counts ($(length(contents))!=$length(counts))."
        return new{contents|>eltype}(form|>uppercase,contents,[0,cumsum(counts)...])
    end
    function AbelianNumbers(form::Char,contents::Vector{<:AbelianNumber},indptr::Vector{Int},::QnsIndptr)
        @assert form|>uppercase ∈ ('G','U','C') "AbelianNumbers error: 'form'($form) is not 'G','U' or 'C'."
        @assert length(contents)+1==length(indptr) "AbelianNumbers error: dismatch shapes of contents and indptr ($(length(contents))+1!=$length(indptr))."
        return new{contents|>eltype}(form|>uppercase,contents,indptr)
    end
end

"""
    AbelianNumbers(qn::AbelianNumber,count::Int=1)

Construct a `AbelianNumbers` with one unique quantum number which occurs `count` times.
"""
AbelianNumbers(qn::AbelianNumber,count::Int=1)=AbelianNumbers('C',[qn],[0,count],qnsindptr)

"""
    AbelianNumbers(od::OrderedDict{<:AbelianNumber,Int})

Construct a `AbelianNumbers` from an ordered dict containing concrete quantum numbers and their counts.
"""
function AbelianNumbers(od::OrderedDict{<:AbelianNumber,Int})
    contents=Vector{od|>keytype}(undef,length(od))
    indptr=zeros(Int,length(od)+1)
    for (i,(qn,count)) in enumerate(od)
        @inbounds contents[i]=qn
        @inbounds indptr[i+1]=indptr[i]+count
    end
    AbelianNumbers('U',contents,indptr,qnsindptr)
end

"""
    AbelianNumbers(od::OrderedDict{<:AbelianNumber,UnitRange{Int}})

Construct a `AbelianNumbers` from an ordered dict containing concrete quantum numbers and their slices.
"""
function AbelianNumbers(od::OrderedDict{<:AbelianNumber,UnitRange{Int}})
    contents=Vector{od|>keytype}(undef,length(od))
    indptr=zeros(Int,length(od)+1)
    for (i,(qn,slice)) in enumerate(od)
        @inbounds contents[i]=qn
        @inbounds @assert indptr[i]+1==slice.start "AbelianNumbers error: slice not consistent."
        @inbounds indptr[i+1]=slice.stop
    end
    AbelianNumbers('U',contents,indptr,qnsindptr)
end

"""
    ==(qns1::AbelianNumbers,qns2::AbelianNumbers) -> Bool
    isequal(qns1::AbelianNumbers,qns2::AbelianNumbers) -> Bool

Overloaded equivalent operator. Two `AbelianNumbers`es are equal to each other if and only if both their `contents`es and `indptr`s are elementwise equal to each other.
!!! note
    It is not necessary for two `AbelianNumbers`es to have the same eltype nor the same form to be equal to each other.
"""
Base.:(==)(qns1::AbelianNumbers,qns2::AbelianNumbers)=qns1.contents==qns2.contents && qns1.indptr==qns2.indptr
Base.isequal(qns1::AbelianNumbers,qns2::AbelianNumbers)=isequal(qns1.contents,qns2.contents) && isequal(qns1.indptr,qns2.indptr)

"""
    show(io::IO,qns::AbelianNumbers)

Show a `AbelianNumbers`.
"""
Base.show(io::IO,qns::AbelianNumbers)=@printf io "QNS(%s)" join((@sprintf("%s=>%s",qn,slice) for (qn,slice) in pairs(qns,qnsindptr)),',')

"""
    string(qns::AbelianNumbers) -> String

Convert a `AbelianNumbers` to string.
"""
Base.string(qns::AbelianNumbers)=@sprintf "QNS(%s,%s)" qns|>length qns|>dimension

"""
    length(qns::AbelianNumbers) -> Int

Get the number of unduplicate qunatum numbers in the `AbelianNumbers`.
"""
Base.length(qns::AbelianNumbers)=length(qns.contents)

"""
    eltype(::Type{<:AbelianNumbers{QN}}) where {QN<:AbelianNumber}
    eltype(qns::AbelianNumbers)

Get the type of the concrete `AbelianNumber` contained in a `AbelianNumbers`.
"""
Base.eltype(::Type{<:AbelianNumbers{QN}}) where {QN<:AbelianNumber}=QN
Base.eltype(qns::AbelianNumbers)=qns|>typeof|>eltype

"""
    getindex(qns::AbelianNumbers,index::Int) -> eltype(qns)
    getindex(qns::AbelianNumbers,slice::UnitRange{Int}) -> AbelianNumbers
    getindex(qns::AbelianNumbers,indices::Vector{Int}) -> AbelianNumbers

Overloaded `[]` operator.
!!! note
    1. For a `AbelianNumbers`, all these `getindex` functions act on its `contents`, i.e. its compressed data, but not on its expansion, i.e. the uncompressed data. This definition is consistent with the [`length`](@ref) function.
    2. When the index is an integer, the result is a `AbelianNumber`, while when the index is a unit range or a vector of intgers, the result is a `AbelianNumbers`. The logic is quite reasonable because such behaviors are much alike to those of a vector container.
"""
Base.getindex(qns::AbelianNumbers,index::Int)=qns.contents[index]
function Base.getindex(qns::AbelianNumbers,slice::UnitRange{Int})
    contents=qns.contents[slice]
    indptr=qns.indptr[slice.start:slice.stop+1]
    indptr.=indptr.-qns.indptr[slice.start]
    return AbelianNumbers(qns.form,contents,indptr,qnsindptr)
end
function Base.getindex(qns::AbelianNumbers,indices::Vector{Int})
    contents=Vector{qns|>eltype}(undef,length(indices))
    indptr=zeros(Int,length(indices)+1)
    for (i,index) in enumerate(indices)
        contents[i]=qns.contents[index]
        indptr[i+1]=indptr[i]+qns.indptr[index+1]-qns.indptr[index]
    end
    return AbelianNumbers('U',contents,indptr,qnsindptr)
end

"""
    iterate(qns::AbelianNumbers,state::Int=1)
    iterate(rv::Iterators.Reverse{<:AbelianNumbers},state::Int=length(rv.itr,false))

Iterate or reversely iterate over the concrete `AbelianNumber`s contained in a `AbelianNumbers`.
"""
Base.iterate(qns::AbelianNumbers,state::Int=1)=state>length(qns) ? nothing : (@inbounds(qns.contents[state]),state+1)
Base.iterate(rv::Reverse{<:AbelianNumbers},state::Int=length(rv.itr))=state<1 ? nothing : (@inbounds(rv.itr.contents[state]),state-1)

"""
    keys(qns::AbelianNumbers) -> Vector{qns|>eltype}

Iterate over the concrete `AbelianNumber`s contained in a `AbelianNumbers`.
"""
Base.keys(qns::AbelianNumbers)=qns.contents

"""
    values(qns::AbelianNumbers,::QnsIndptr)
    values(qns::AbelianNumbers,::QnsCounts)

Iterate over the slices/counts of the `AbelianNumbers`.
"""
@views Base.values(qns::AbelianNumbers,::QnsIndptr)=((start+1):stop for (start,stop) in zip(qns.indptr[1:end-1],qns.indptr[2:end]))
@views Base.values(qns::AbelianNumbers,::QnsCounts)=(stop-start for (start,stop) in zip(qns.indptr[1:end-1],qns.indptr[2:end]))

"""
    pairs(qns::AbelianNumbers,choice::Union{QnsIndptr,QnsCounts})

Iterate over the `AbelianNumber=>slice` or `AbelianNumber=>count` pairs.
"""
Base.pairs(qns::AbelianNumbers,choice::Union{QnsIndptr,QnsCounts})=Base.Generator(=>,keys(qns),values(qns,choice))

"""
    toordereddict(qns::AbelianNumbers,::QnsIndptr) -> OrderedDict{qns|>eltype,UnitRange{Int}}
    toordereddict(qns::AbelianNumbers,::QnsCounts) -> OrderedDict{qns|>eltype,Int}

Convert a `AbelianNumbers` to an ordered dict.
"""
function toordereddict(qns::AbelianNumbers,::QnsIndptr)
    @assert qns.form != 'G' "toordereddict error: input `AbelianNumbers` cannot be `G` formed."
    result=OrderedDict{qns|>eltype,UnitRange{Int}}()
    for i=1:length(qns)
        @inbounds result[qns.contents[i]]=qns.indptr[i]+1:qns.indptr[i+1]
    end
    return result
end
function toordereddict(qns::AbelianNumbers,::QnsCounts)
    @assert qns.form != 'G' "toordereddict error: input `AbelianNumbers` cannot be `G` formed."
    result=OrderedDict{qns|>eltype,Int}()
    for i=1:length(qns)
        @inbounds result[qns.contents[i]]=qns.indptr[i+1]-qns.indptr[i]
    end
    return result
end

"""
    dimension(qns::AbelianNumbers) -> Int

The dimension of the Hilbert space a `AbelianNumbers` represents.
"""
dimension(qns::AbelianNumbers)=@inbounds(qns.indptr[end])

"""
    sort(qns::AbelianNumbers) -> AbelianNumbers,Vector{Int}

Sort the quantum numbers of a `AbelianNumber`, return the sorted `AbelianNumber` and the permutation array that sorts the expansion of the original `AbelianNumbers`.
"""
function Base.sort(qns::AbelianNumbers)
    ctpts=sortperm(qns.contents,alg=Base.Sort.QuickSort)
    masks=Vector{Bool}(undef,length(qns.contents));masks[1]=true
    unduplicate=1
    for i=2:length(ctpts)
        @inbounds masks[i]=qns.contents[ctpts[i]]≠qns.contents[ctpts[i-1]]
        @inbounds masks[i] && (unduplicate+=1)
    end
    contents=Vector{qns|>eltype}(undef,unduplicate)
    indptr=zeros(Int,unduplicate+1)
    permutation=Vector{Int}(undef,dimension(qns))
    qncount,ptcount=0,0
    for (mask,index) in zip(masks,ctpts)
        @inbounds mask && (qncount+=1;contents[qncount]=qns.contents[index];indptr[qncount+1]=indptr[qncount])
        @inbounds indptr[qncount+1]+=qns.indptr[index+1]-qns.indptr[index]
        for (i,p) in enumerate(@inbounds(qns.indptr[index]+1:qns.indptr[index+1]))
            @inbounds permutation[ptcount+i]=p
        end
        @inbounds ptcount+=qns.indptr[index+1]-qns.indptr[index]
    end
    return AbelianNumbers('C',contents,indptr,qnsindptr),permutation
end

struct QnsCompression   <: QnsProtocol end
struct QnsExpansion     <: QnsProtocol end
"""
    qnscompression

Indicate that [`findall`](@ref) and [`permute`](@ref) use the compressed contents.
"""
const qnscompression=QnsCompression()
"""
    qnsexpansion

Indicate that [`findall`](@ref) and [`permute`](@ref) use the expanded contents.
"""
const qnsexpansion=QnsExpansion()

"""
    findall(target::QN,qns::AbelianNumbers{QN},::QnsCompression) where QN<:AbelianNumber -> Vector{Int}
    findall(target::QN,qns::AbelianNumbers{QN},::QnsExpansion) where QN<:AbelianNumber -> Vector{Int}
    findall(targets::NTuple{N,QN},qns::AbelianNumbers{QN},::QnsCompression) where {N,QN<:AbelianNumber} -> Vector{Int}
    findall(targets::NTuple{N,QN},qns::AbelianNumbers{QN},::QnsExpansion) where {N,QN<:AbelianNumber} -> Vector{Int}

Find all the indices of the target quantum numbers in the contents ([`qnscompression`](@ref) case) or the expansion ([`qnsexpansion`](@ref) case) of a `AbelianNumbers`.
"""
Base.findall(target::QN,qns::AbelianNumbers{QN},::QnsCompression) where QN<:AbelianNumber=findall((target,),qns,QnsCompression())
Base.findall(target::QN,qns::AbelianNumbers{QN},::QnsExpansion) where QN<:AbelianNumber=findall((target,),qns,QnsExpansion())
function Base.findall(targets::NTuple{N,QN},qns::AbelianNumbers{QN},::QnsCompression) where {N,QN<:AbelianNumber}
    result=Int[]
    if qns.form=='C'
        for qn in targets
            range=searchsorted(qns.contents,qn)
            range.start<=range.stop && push!(result,range.start)
        end
    else
        for (i,qn) in enumerate(qns)
            qn ∈ targets && push!(result,i)
        end
    end
    result
end
function Base.findall(targets::NTuple{N,QN},qns::AbelianNumbers{QN},::QnsExpansion) where {N,QN<:AbelianNumber}
    result=Int[]
    for index in findall(targets,qns,qnscompression)
        @inbounds append!(result,(qns.indptr[index]+1):qns.indptr[index+1])
    end
    result
end

"""
    filter(target::QN,qns::AbelianNumbers{QN}) where QN<:AbelianNumber -> AbelianNumbers{QN}
    filter(targets::NTuple{N,QN},qns::AbelianNumbers{QN}) where {N,QN<:AbelianNumber} -> AbelianNumbers{QN}

Find a subset of a `AbelianNumbers` by picking out the quantum numbers in targets.
"""
Base.filter(target::QN,qns::AbelianNumbers{QN}) where QN<:AbelianNumber=qns[findall(target,qns,qnscompression)]
Base.filter(targets::NTuple{N,QN},qns::AbelianNumbers{QN}) where {N,QN<:AbelianNumber}=qns[findall(targets,qns,qnscompression)]

"""
    permute(qns::AbelianNumbers,permutation::Vector{Int},::QnsCompression) -> AbelianNumbers
    permute(qns::AbelianNumbers,permutation::Vector{Int},::QnsExpansion) -> AbelianNumbers

Reorder the quantum numbers contained in a `AbelianNumbers` with a permutation and return the new one.

For [`qnscompression`](@ref) case, the permutation is for the compressed contents of the original `AbelianNumbers` while for [`qnsexpansion`](@ref) case, the permutation is for the expanded contents of the original `AbelianNumbers`.
"""
permute(qns::AbelianNumbers,permutation::Vector{Int},::QnsCompression)=qns[permutation]
function permute(qns::AbelianNumbers,permutation::Vector{Int},::QnsExpansion)
    contents=Vector{qns|>eltype}(undef,length(permutation))
    indptr=zeros(Int,length(permutation)+1)
    expansion=expand(qns,qnsindices)
    for (i,p) in enumerate(permutation)
        contents[i]=qns.contents[expansion[p]]
        indptr[i+1]=indptr[i]+1
    end
    return AbelianNumbers('G',contents,indptr,qnsindptr)
end

struct QnsContents  <: QnsProtocol end
struct QnsIndices   <: QnsProtocol end
"""
    qnscontents

Indicate that [`expand`](@ref) uses the compressed/expanded contents.
"""
const qnscontents=QnsContents()
"""
    qnsindices

Indicate that [`expand`](@ref) uses the indices of the compressed/expanded contents.
"""
const qnsindices=QnsIndices()

"""
    expand(qns::AbelianNumbers,::QnsContents) -> Vector{qns|>eltype}
    expand(qns::AbelianNumbers,::QnsIndices) -> Vector{Int}

Expand the contents ([`qnscontents`](@ref) case) or indices ([`qnsindices`](@ref) case) of a `AbelianNumbers` to the uncompressed form.
"""
function expand(qns::AbelianNumbers,::QnsContents)
    result=Vector{qns|>eltype}(undef,dimension(qns))
    for i=1:length(qns)
        for j=@inbounds((qns.indptr[i]+1):qns.indptr[i+1])
            @inbounds result[j]=qns.contents[i]
        end
    end
    return result
end
function expand(qns::AbelianNumbers,::QnsIndices)
    result=Vector{Int}(undef,dimension(qns))
    for i=1:length(qns)
        for j=@inbounds((qns.indptr[i]+1):qns.indptr[i+1])
            @inbounds result[j]=i
        end
    end
    return result
end

"""
    +(qn::AbelianNumber) -> typeof(qn)
    +(qn::QN,qns::QN...) where QN<:AbelianNumber -> QN
    +(qns::AbelianNumbers) -> AbelianNumbers
    +(qn::QN,qns::AbelianNumbers{QN}) where QN<:AbelianNumber -> AbelianNumbers{QN}
    +(qns::AbelianNumbers{QN},qn::QN) where QN<:AbelianNumber -> AbelianNumbers{QN}

Overloaded `+` operator for `AbelianNumber` and `AbelianNumbers`.
!!! note
    1. The addition between a `AbelianNumbers` and an `AbelianNumber` is just a global shift of the contents of the `AbelianNumbers` by the `AbelianNumber`, therefore, the result is a `AbelianNumbers`.
    2. `+` cannot be used between two `AbelianNumbers` because the result is ambiguous. Instead, use `⊕` for direct sum and `⊗` for direct product.
    3. To ensure type stability, two `AbelianNumber` can be added together if and only if they are of the same type.
    4. Similarly, a `AbelianNumber` and a `AbelianNumbers` can be added together if and only if the former's type is the same with the latter's eltype.
"""
Base.:+(qn::AbelianNumber)=qn
Base.:+(qn::QN,qns::QN...) where QN<:AbelianNumber=map(+,qn,qns...)
Base.:+(qns::AbelianNumbers)=qns
Base.:+(qn::QN,qns::AbelianNumbers{QN}) where QN<:AbelianNumber=qns+qn
Base.:+(qns::AbelianNumbers{QN},qn::QN) where QN<:AbelianNumber=AbelianNumbers(qns.form=='G' ? 'G' : 'U',[iqn+qn for iqn in qns],qns.indptr,qnsindptr)

"""
    -(qn::AbelianNumber) -> typeof(qn)
    -(qn1::QN,qn2::QN) where QN<:AbelianNumber -> QN
    -(qns::AbelianNumbers) -> AbelianNumbers
    -(qn::QN,qns::AbelianNumbers{QN}) where QN<:AbelianNumber -> AbelianNumbers{QN}
    -(qns::AbelianNumbers{QN},qn::QN) where QN<:AbelianNumber -> AbelianNumbers{QN}

Overloaded `-` operator for `AbelianNumber` and `AbelianNumbers`.
!!! note
    1. The subtraction between a `AbelianNumbers` and a `AbelianNumber` is just a global shift of the contents of the `AbelianNumbers` by the `AbelianNumber`, therefore, the result is a `AbelianNumbers`.
    2. `-` cannot be used between two `AbelianNumbers` because the result is ambiguous. Instead, use `⊕` with signs for direct sum and `⊗` with signs for direct product.
    3. To ensure type stability, a `AbelianNumber` can be subtracted by another `AbelianNumber` if and only if they are of the same type.
    4. Similarly, a `AbelianNumber` can be subtracted by a `AbelianNumbers` or vice versa if and only if the former's type is the same with the latter's eltype.
"""
Base.:-(qn::AbelianNumber)=map(-,qn)
Base.:-(qn1::QN,qn2::QN) where QN<:AbelianNumber=map(-,qn1,qn2)
Base.:-(qns::AbelianNumbers)=AbelianNumbers(qns.form=='G' ? 'G' : 'U',-qns.contents,qns.indptr,qnsindptr)
Base.:-(qn::QN,qns::AbelianNumbers{QN}) where QN<:AbelianNumber=AbelianNumbers(qns.form=='G' ? 'G' : 'U',[qn-iqn for iqn in qns],qns.indptr,qnsindptr)
Base.:-(qns::AbelianNumbers{QN},qn::QN) where QN<:AbelianNumber=AbelianNumbers(qns.form=='G' ? 'G' : 'U',[iqn-qn for iqn in qns],qns.indptr,qnsindptr)

"""
    *(qn::AbelianNumber,factor::Integer) -> typeof(qn)
    *(factor::Integer,qn::AbelianNumber) -> typeof(qn)
    *(qns::AbelianNumbers,factor::Integer) -> AbelianNumbers
    *(factor::Integer,qns::AbelianNumbers) -> AbelianNumbers

Overloaded `*` operator for the multiplication between an integer and a `AbelianNumber` or a `AbelianNumbers`.
"""
@generated function Base.:*(qn::AbelianNumber,factor::Integer)
    exprs=Expr[:(getfield(qn,$i)*factor) for i=1:(qn|>fieldnames|>length)]
    return :(typeof(qn)($(exprs...)))
end
Base.:*(factor::Integer,qn::AbelianNumber)=qn*factor
Base.:*(qns::AbelianNumbers,factor::Integer)=AbelianNumbers(qns.form=='G' ? 'G' : 'U',[qn*factor for qn in qns.contents],qns.indptr,qnsindptr)
Base.:*(factor::Integer,qns::AbelianNumbers)=qns*factor

"""
    ^(qn::AbelianNumber,factor::Integer) -> typeof(qn)
    ^(qns::AbelianNumbers,factor::Integer) -> AbelianNumbers

Overloaded `^` operator for `AbelianNumber` and `AbelianNumbers`. This operation translates into the direct product of `factor` copies of `qn` or `qns`.
"""
Base.:^(qn::AbelianNumber,factor::Integer)=kron(NTuple{factor,qn|>typeof}(qn for i=1:factor)...)
Base.:^(qns::AbelianNumbers,factor::Integer)=kron(NTuple{factor,qns|>typeof}(qns for i=1:factor)...)

"""
    ⊕(qns::AbelianNumber...) -> AbelianNumbers{qns|>eltype}
    ⊕(qnses::AbelianNumbers...) -> qnses|>eltype

Get the direct sum of some `AbelianNumber`s or `AbelianNumbers`es.
"""
⊕(qns::AbelianNumber...)=union(qns...)
⊕(qnses::AbelianNumbers...)=union(qnses...)

"""
    ⊗(qns::AbelianNumber...) -> eltype(qns)
    ⊗(qnses::AbelianNumbers...) -> eltype(qnses)

Get the direct product of some `AbelianNumber`s or `AbelianNumbers`es.
"""
⊗(qns::AbelianNumber...)=kron(qns...)
⊗(qnses::AbelianNumbers...)=kron(qnses...)

"""
    union(qns::Vararg{<:AbelianNumber,N};signs::NTuple{N,Int}=ntuple(i->1,N)) where N -> AbelianNumbers
    union(qnses::Vararg{AbelianNumbers{QN},N};signs::NTuple{N,Int}=ntuple(i->1,N)) where {N,QN<:AbelianNumber} -> AbelianNumbers{QN}

Get the direct sum of some `AbelianNumber`s or `AbelianNumbers`es.
!!! note
    1. Physically, the direct sum of a couple of `AbelianNumber`s or `AbelianNumbers`es is defined by the direct sum of the bases of the Hilbert spaces they represent. Therefore, the input `AbelianNumber`s or `AbelianNumbers`es must be homogenous. Inhomogenous 'AbelianNumber's must be direct producted first to ensure homogenity before the direct sum.
    2. Apparently, the dimension of the result equals the summation of those of the inputs, which means, even for `AbelianNumber`s, the result will be naturally a `AbelianNumbers` because the dimension of the result is larger than 1.
    3. Signs of `AbelianNumber`s or `AbelianNumbers`es can be provided when getting their direct sums.
"""
function Base.union(qns::Vararg{<:AbelianNumber,N};signs::NTuple{N,Int}=ntuple(i->1,N)) where N
    AbelianNumbers('G',[sign==1 ? qn : -qn for (sign,qn) in zip(signs,qns)],fill(1,qns|>length),qnscounts)
end
function Base.union(qnses::Vararg{AbelianNumbers{QN},N};signs::NTuple{N,Int}=ntuple(i->1,N)) where {N,QN<:AbelianNumber}
    lengths=NTuple{N,Int}(length(qns) for qns in qnses)
    contents=Vector{QN}(undef,sum(lengths))
    indptr=zeros(Int,sum(lengths)+1)
    count=0
    for (sign,qns,length) in zip(signs,qnses,lengths)
        for i=1:length
            @inbounds contents[count+i]=sign==1 ? qns.contents[i] : -qns.contents[i]
            @inbounds indptr[count+i+1]=indptr[count+i]+qns.indptr[i+1]-qns.indptr[i]
        end
        count+=length
    end
    AbelianNumbers('G',contents,indptr,qnsindptr)
end

"""
    kron(::Type{QN},qn1::AbelianNumber,qn2::AbelianNumber) where QN<:AbelianNumber -> QN
    kron(qns::Vararg{<:AbelianNumber,N};signs::NTuple{N,Int}=ntuple(i->1,N)) where N -> eltype(qns)
    kron(qnses::Vararg{AbelianNumbers{QN},N};signs::NTuple{N,Int}=ntuple(i->1,N)) where {N,QN<:AbelianNumber} -> AbelianNumbers{QN}

Get the direct product of some `AbelianNumber`s or `AbelianNumbers`es.
!!! note
    1. Physically, the direct product of a couple of `AbelianNumber`s or `AbelianNumbers`es are defined by the direct product of the bases of the Hilbert spaces they represent. Therefore, `AbelianNumbers` with differenct types or `AbelianNumbers`es with differenct eltypes are allowed to be direct producted in principle. However, for simplicity, we only implement a method which handle the situation of two `AbelianNumber`s with differenct types. The type of the result should be provided as the first parameter. Note that in this situation, the `fieldnames` and `periods` of the result type must be exactly equal to the flattened fieldnames and periods of the two input `AbelianNumber`s, which means, even the order of the input `AbelianNumber`s matters.
    2. Apparently, the dimension of the result equals the product of those of the inputs. Therefore, the direct product of `AbelianNumber`s is also a `AbelianNumber` since its dimension is still one.
    3. For other situations except the one mentioned in Note.1, the input `AbelianNumber`s or `AbelianNumbers`es must be homogenous. Meanwhile, signs can also be provided for these situations. Note that each quantum number in the contents of the result is obtained by a summation of the corresponding quanum numbers out of the inputs with the correct signs. This is a direct observation of the Abelian nature of our quantum numbers.
"""
function Base.kron(::Type{QN},qn1::AbelianNumber,qn2::AbelianNumber) where QN<:AbelianNumber
    qnnames,qn1names,qn2names=QN|>fieldnames,qn1|>typeof|>fieldnames,qn2|>typeof|>fieldnames
    qnperiods,qn1periods,qn2periods=QN|>periods,qn1|>typeof|>periods,qn2|>typeof|>periods
    @assert qnnames==NTuple{QN|>length,Symbol}(flatten((qn1names,qn2names))) "kron error: fieldnames not match ($(qnnames),$(qn1names),$(qn2names))."
    @assert qnperiods==NTuple{QN|>length,Float}(flatten((qn1periods,qn2periods))) "kron error: periods not match ($(qnperiods),$(qn1periods),$(qn2periods))."
    QN(qn1...,qn2...)
end
Base.kron(qns::Vararg{<:AbelianNumber,N};signs::NTuple{N,Int}=ntuple(i->1,N)) where N=sum(sign==1 ? qn : -qn for (sign,qn) in zip(signs,qns))
function Base.kron(qnses::Vararg{AbelianNumbers{QN},N};signs::NTuple{N,Int}=ntuple(i->1,N)) where {N,QN<:AbelianNumber}
    lengths=NTuple{N,Int}(i<N ? dimension(qns) : length(qns) for (i,qns) in enumerate(qnses))
    contents=Vector{QN}(undef,prod(lengths))
    indptr=zeros(Int,prod(lengths)+1)
    cache=Vector{QN}(undef,N)
    @inbounds expansions=NTuple{N-1,Vector{Int}}(expand(qnses[i],qnsindices) for i=1:(N-1))
    for (i,indices) in enumerate(product(NTuple{N,UnitRange{Int64}}(1:length for length in reverse(lengths))...))
        for (j,(sign,qns,index)) in enumerate(zip(signs,qnses,reverse(indices)))
            @inbounds pos=j<N ? expansions[j][index] : index
            @inbounds cache[j]=sign==1 ? qns.contents[pos] : -qns.contents[pos]
        end
        @inbounds indptr[i+1]=indptr[i]+qnses[end].indptr[indices[1]+1]-qnses[end].indptr[indices[1]]
        @inbounds contents[i]=sum(cache)
    end
    AbelianNumbers('G',contents,indptr,qnsindptr)
end

"""
    ukron(qnses::Vararg{AbelianNumbers{QN},N};signs::NTuple{N,Int}=ntuple(i->1,N)) where {N,QN<:AbelianNumber} -> AbelianNumbers{QN},Dict{QN,Dict{NTuple{N,QN},UnitRange{Int}}}

Unitary Kronecker product of several `AbelianNumbers`es. The product result as well as the records of the product will be returned.
!!! note
    1. All input `AbelianNumbers` must be 'U' formed or 'C' formed.
    2. Since duplicate quantum number are not allowed in 'U' formed and 'C' formed `AbelianNumbers`es, in general, there exists a merge process of duplicate quantum numbers in the product result. Therefore, records are needed to keep track of this process, which will be returned along with the product result. The records are stored in a `Dict{QN,Dict{NTuple{N,QN},UnitRange{Int}}}` typed dict, in which, for each unduplicate quantum number `qn` in the product result, there exist a record `Dict((qn₁,qn₂,...)=>start:stop,...)` telling what quantum numbers `(qn₁,qn₂,...)` a mereged duplicate `qn` comes from and what slice `start:stop` this merged duplicate corresponds.
"""
function ukron(qnses::Vararg{AbelianNumbers{QN},N};signs::NTuple{N,Int}=ntuple(i->1,N)) where {N,QN<:AbelianNumber}
    @assert all(qns.form=='U' || qns.form=='C' for qns in qnses) "ukron error: all input qnses should be 'U' formed or 'C' formed."
    lengths=NTuple{N,Int}(length(qns) for qns in qnses)
    cache=Vector{QN}(undef,N)
    container=OrderedDict{QN,Int}()
    records=Dict{QN,Dict{NTuple{N,QN},UnitRange{Int}}}()
    for indices in product(NTuple{N,UnitRange{Int64}}(1:length for length in reverse(lengths))...)
        qn,count=QN|>zero,1
        for (j,(sign,qns,index)) in enumerate(zip(signs,qnses,reverse(indices)))
            @inbounds cache[j]=qns.contents[index]
            @inbounds qn=(sign==1 ? qn+cache[j] : qn-cache[j])
            @inbounds count=count*(qns.indptr[index+1]-qns.indptr[index])
        end
        container[qn]=get(container,qn,0)+count
        !haskey(records,qn) && (records[qn]=Dict{NTuple{N,QN},UnitRange{Int}}())
        records[qn][NTuple{N,QN}(cache)]=(container[qn]-count+1):container[qn]
    end
    contents=QN[qn for qn in keys(container)]
    sort!(contents)
    indptr=zeros(Int,length(container)+1)
    for (i,qn) in enumerate(contents)
        @inbounds indptr[i+1]=indptr[i]+container[qn]
    end
    AbelianNumbers('C',contents,indptr,qnsindptr),records
end

struct QnsBruteForce <: QnsProtocol end
struct QnsMonteCarlo <: QnsProtocol end
"""
    qnsbruteforce

Indicate that [`decompose`](@ref) uses the brute force method.
"""
const qnsbruteforce=QnsBruteForce()

"""
    qnsmontecarlo

Indicate that [`decompose`](@ref) uses the Monte Carlo method.
"""
const qnsmontecarlo=QnsMonteCarlo()

"""
    decompose(qnses::NTuple{N,AbelianNumbers{QN}},target::QN,signs::NTuple{N,Int},::QnsBruteForce;nmax::Int=20) where {N,QN<:AbelianNumber} -> Vector{NTuple{N,Int}}
    decompose(qnses::NTuple{N,AbelianNumbers{QN}},target::QN,signs::NTuple{N,Int},::QnsMonteCarlo;nmax::Int=20) where {N,QN<:AbelianNumber} -> Vector{NTuple{N,Int}}

Find a couple of decompositions of `target` with respect to `qnses`.
!!! note
    A tuple of integers `(i₁,i₂,...)` is called a decomposition of a given `target` with respect to the given `qnses` if and only if they satisfy the "decomposition rule":
    ```math
    \\sum_\\text{j} \\text{signs}[\\text{j}]\\times\\text{qnses}[\\text{j}][\\text{i}_{\\text{j}}]==\\text{target}
    ```
    This equation is in fact a kind of a set of restricted [linear Diophantine equations](https://en.wikipedia.org/wiki/Diophantine_equation#Linear_Diophantine_equations). Indeed, our quantum numbers are always discrete Abelian ones and all instances of a concrete `AbelianNumber` forms a [module](https://en.wikipedia.org/wiki/Module_(mathematics)) over the [ring](https://en.wikipedia.org/wiki/Ring_(mathematics)) of integers. Therefore, each quantum number can be represented as a integral multiple of the unit element of the Abelian module, which results in the final reduction of the above equation to a set of linear Diophantine equations. Then finding a decomposition is equivalent to find a solution of the reduced linear Diophantine equations, with the restriction that the quantum numbers constructed from the solution should be in the corresponding `qnses`. Here we provide two methods to find such decompositions, one is by brute force ([`qnsbruteforce`](@ref) case), and the other is by Monte Carlo simultatioins ([`qnsmontecarlo`](@ref) case).
"""
function decompose(qnses::NTuple{N,AbelianNumbers{QN}},target::QN,signs::NTuple{N,Int},::QnsBruteForce;nmax::Int=20) where {N,QN<:AbelianNumber}
    result=Set{NTuple{N,Int}}()
    cache=Vector{Int}(undef,N)
    dimensions=NTuple{N,Int}(dimension(qns) for qns in reverse(qnses))
    indices=findall(target,kron(qnses...,signs=signs),qnsexpansion)
    nmax<length(indices) && (shuffle!(MersenneTwister(),indices);indices=@views indices[1:nmax])
    for index in indices
        for (i,dimension) in enumerate(dimensions)
            @inbounds cache[end+1-i]=(index-1)%dimension+1
            @inbounds index=(index-1)÷dimension+1
        end
        push!(result,NTuple{N,Int}(cache))
    end
    return collect(NTuple{N,Int},result)
end
function decompose(qnses::NTuple{N,AbelianNumbers{QN}},target::QN,signs::NTuple{N,Int},::QnsMonteCarlo;nmax::Int=20) where {N,QN<:AbelianNumber}
    seed!()
    result=Set{NTuple{N,Int}}()
    expansions=Vector{QN}[expand(qns,qnscontents) for qns in qnses]
    dimensions=NTuple{N,Int}(dimension(qns) for qns in qnses)
    diff=indices->norm(sum(sign==1 ? expansion[index] : -expansion[index] for (sign,expansion,index) in zip(signs,expansions,indices))-target)
    count=1
    while true
        oldindices,newindices=Int[rand(1:dimension) for dimension in dimensions],ones(Int,N)
        olddiff,newdiff=diff(oldindices),diff(newindices)
        while newdiff>0
            pos=rand(1:N)
            index=rand(1:dimensions[pos]-1)
            @inbounds newindices[pos]=index<newindices[pos] ? index : index+1
            newdiff=diff(newindices)
            if newdiff<=olddiff || exp(olddiff-newdiff)>rand()
                oldindices[:]=newindices
                olddiff=newdiff
            end
            newindices[:]=oldindices
        end
        count=count+1
        push!(result,NTuple{N,Int}(newindices))
        (length(result)>=nmax || count>nmax*5) && break
    end
    return collect(NTuple{N,Int},result)
end

"""
    SQNS(S::Real) -> AbelianNumbers{SQN}

Construct the `AbelianNumbers` of the Hilbert space of a signle spin `S`.
"""
SQNS(S::Real)=AbelianNumbers('C',[SQN(sz) for sz=-S:S],collect(0:Int(2*S+1)),qnsindptr)

"""
    PQNS(N::Real) -> AbelianNumbers{PQN}

Construct the `AbelianNumbers` of the Hilbert space of a single-particle state with at most `N` identical particles.
"""
PQNS(N::Real)=AbelianNumbers('C',[PQN(np) for np=0:N],collect(0:Int(N)+1),qnsindptr)

"""
    SzPQNS(Sz::Real) -> AbelianNumbers{SPQN}

Construct the `AbelianNumbers` of the Hilbert space of a single-paritcle state with at most one particle whose spin-z component is `Sz`.
"""
SzPQNS(Sz::Real)=AbelianNumbers('C',[SPQN(0.0,0.0),SPQN(1.0,Sz)],[0,1,2],qnsindptr)

"""
    SPQNS(S::Real) -> AbelianNumbers{SPQN}

Construct the `AbelianNumbers` of the Hilbert space of a single site with internal degrees of freedom that can be ascribed to a spin `S`.
"""
function SPQNS(S::Real)
    sqns=[SQN(sz) for sz=-S:S]
    contents=[SPQN(0.0,0.0)]
    for n=1:length(sqns)
        pn=PQN(n*1.0)
        for spins in Combinations{n}(sqns)
            push!(contents,kron(SPQN,pn,sum(spins)))
        end
    end
    return sort(AbelianNumbers('G',contents,collect(0:length(contents)),qnsindptr))[1]
end

"""
    Z2QNS() -> AbelianNumbers{Z2QN}

Construct the `AbelianNumbers` of a ``Z_2`` Hilbert space.
"""
Z2QNS()=AbelianNumbers('C',[Z2QN(0.0),Z2QN(1.0)],[0,1,2],qnsindptr)

end #module
