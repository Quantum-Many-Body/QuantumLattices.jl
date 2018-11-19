module QuantumNumber

using Base.Iterators: Reverse,flatten,product,reverse
using Printf: @printf,@sprintf
using DataStructures: OrderedDict
using Random: MersenneTwister,seed!,shuffle!
using LinearAlgebra: norm
using Combinatorics: combinations
using ..NamedVector: HomoNamedVector

export AbstractQuantumNumber
export regularize!,regularize
export @quantumnumber,periods,SQN,PQN,SPQN,Z2QN
export qnscounts,qnsindptr,qnscontents,qnsexpansion,qnsindices,qnsbruteforce,qnsmontecarlo
export QuantumNumbers
export ⊕,⊗,ukron,dimension,expand,decompose,subset,reorder,toordereddict
export SQNS,PQNS,SzPQNS,SPQNS,Z2QNS

"Abstract type for all concrete quantum numbers for a single basis."
abstract type AbstractQuantumNumber <: HomoNamedVector{Float64} end

"""
    dimension(::Type{<:AbstractQuantumNumber}) -> Int
    dimension(::AbstractQuantumNumber) -> Int

The dimension of the Hilbert space a `AbstractQuantumNumber` represents. Apparently, this is always 1.
"""
dimension(::Type{<:AbstractQuantumNumber})=1
dimension(::AbstractQuantumNumber)=1

"""
    regularize!(::Type{QN},array::AbstractVector{Float64}) where QN<:AbstractQuantumNumber
    regularize!(::Type{QN},array::AbstractMatrix{Float64}) where QN<:AbstractQuantumNumber

Regularize the elements of an array in place so that it can represent quantum numbers.
"""
function regularize!(::Type{QN},array::AbstractVector{Float64}) where QN<:AbstractQuantumNumber
    @assert size(array,1)==QN|>length "regularize! error: not consistent shape of input array and $QN."
    for (i,period) in enumerate(QN|>periods)
        @inbounds period!=Inf && (array[i]=array[i]%period+(array[i]%period<0 ? period : 0))
    end
end
function regularize!(::Type{QN},array::AbstractMatrix{Float64}) where QN<:AbstractQuantumNumber
    @assert size(array,1)==QN|>length "regularize! error: not consistent shape of input array and $QN."
    for (i,period) in enumerate(QN|>periods)
        if period!=Inf
            for j in 1:size(array,2)
                @inbounds array[i,j]=array[i,j]%period+(array[i,j]%period<0 ? period : 0)
            end
        end
    end
end

"""
    regularize(::Type{QN},array::Union{AbstractVector{Float64},AbstractMatrix{Float64}}) where {QN<:AbstractQuantumNumber}

Regularize the elements of an array and return a copy that can represent quantum numbers.
"""
function regularize(::Type{QN},array::Union{AbstractVector{Float64},AbstractMatrix{Float64}}) where {QN<:AbstractQuantumNumber}
    result=copy(array)
    regularize!(QN,result)
    return result
end

"""
    @quantumnumber typename fieldnames fieldperiods

Construct a concrete `AbstractQuantumNumber` with the type name being `typename`, fieldnames specified by `fieldnames` and periods specified by `fieldperiods`.
"""
macro quantumnumber(typename,fieldnames,fieldperiods)
    typename=Symbol(typename)
    fieldnames=tuple(eval(fieldnames)...)
    fieldperiods=tuple(eval(fieldperiods)...)
    arguments=ntuple(i->Symbol(:v,i),length(fieldnames))
    @assert length(fieldnames)==length(fieldperiods) "quantumnumber error: number of fieldnames($(length(fieldnames))) and fieldperiods($(length(fieldperiods))) not equal."
    @assert all(isa(name,Symbol) for name in fieldnames) "quantumnumber error: all field names should be Symbol."
    @assert all(fieldperiods.>0) "quantumnumber error: all field periods should be positive."
    if all(fieldperiods.==Inf)
        title=Expr(:call,:($(esc(typename))),(:($arg::Float64) for arg in arguments)...)
        body=Expr(:call,:new,arguments...)
    else
        title=Expr(:call,:($(esc(typename))),(:($arg::Float64) for arg in arguments)...,Expr(:kw,:(regularize::Bool),:true))
        body=Expr(:call,:new,(p==Inf ? arg : :(regularize ? ($arg)%($p)+(($arg)%($p)<0 ? $p : 0) : $arg) for (p,arg) in zip(fieldperiods,arguments))...)
    end
    newtype=Expr(:struct,false,:($(esc(typename))<:AbstractQuantumNumber),Expr(:block,(:($field::Float64) for field in fieldnames)...,Expr(:(=),title,body)))
    functions=Expr(:block,:(Base.fieldnames(::Type{<:$(esc(typename))})=$fieldnames),:(periods(::Type{<:$(esc(typename))})=$fieldperiods))
    return Expr(:block,:(global periods),:(Base.@__doc__($newtype)),functions)
end

"""
    SQN(Sz::Real)

The concrete `AbstractQuantumNumber` of a quantum system with spin z-component `Sz` conserved.
"""
@quantumnumber "SQN" (:Sz,) (Inf,)

"""
    PQN(N::Real)

The concrete `AbstractQuantumNumber` of a quantum system with particle number `N` conserved.
"""
@quantumnumber "PQN" (:N,) (Inf,)

"""
    SPQN(N::Real,Sz::Real)

The concrete `AbstractQuantumNumber` of a quantum system with both particle number `N` and spin z-component `Sz` conserved.
"""
@quantumnumber "SPQN" (:N,:Sz) (Inf,Inf)

"""
    Z2QN(N::Real)

The concrete `AbstractQuantumNumber` of a quantum system with a Z₂-like conserved quantity.
"""
@quantumnumber "Z2QN" (:N,) (2,)

"Choice associated with `QuantumNumbers`, meaning 'by indptr'."
const qnsindptr=Val(1)

"Choice associated with `QuantumNumbers`, meaning 'by counts'."
const qnscounts=Val(2)

"Choice associated with `QuantumNumbers`, meaning 'for contents'."
const qnscontents=Val(3)

"Choice associated with `QuantumNumbers`, meaning 'for expansion'."
const qnsexpansion=Val(4)

"Choice associated with `QuantumNumbers`, meaning 'for indices'."
const qnsindices=Val(5)

"Choice associated with `QuantumNumbers`, meaning 'by brute force'."
const qnsbruteforce=Val(6)

"Choice associated with `quantumnumbers`, meaning 'by Monte Carlo'."
const qnsmontecarlo=Val(7)

"""
    QuantumNumbers(form::Char,contents::Vector{<:AbstractQuantumNumber},counts::Vector{Int},::typeof(qnscounts))
    QuantumNumbers(form::Char,contents::Vector{<:AbstractQuantumNumber},indptr::Vector{Int},::typeof(qnsindptr))

The whole quantum numbers of the total bases of a Hilbert space. The default constructors construct a `QuantumNumbers` from a vector of concrete quantum numbers and an vector containing their counts or indptr.
"""
struct QuantumNumbers{QN<:AbstractQuantumNumber}
    form::Char
    contents::Vector{QN}
    indptr::Vector{Int}
    function QuantumNumbers(form::Char,contents::Vector{<:AbstractQuantumNumber},counts::Vector{Int},::typeof(qnscounts))
        @assert form|>uppercase ∈ ('G','U','C') "QuantumNumbers error: 'form'($form) is not 'G','U' or 'C'."
        @assert length(contents)==length(counts) "QuantumNumbers error: dismatch lengths of contents and counts ($(length(contents))!=$length(counts))."
        return new{contents|>eltype}(form|>uppercase,contents,[0,cumsum(counts)...])
    end
    function QuantumNumbers(form::Char,contents::Vector{<:AbstractQuantumNumber},indptr::Vector{Int},::typeof(qnsindptr))
        @assert form|>uppercase ∈ ('G','U','C') "QuantumNumbers error: 'form'($form) is not 'G','U' or 'C'."
        @assert length(contents)+1==length(indptr) "QuantumNumbers error: dismatch shapes of contents and indptr ($(length(contents))+1!=$length(indptr))."
        return new{contents|>eltype}(form|>uppercase,contents,indptr)
    end
end

"""
    QuantumNumbers(qn::AbstractQuantumNumber,count::Int=1)

Construct a `QuantumNumbers` with one unique quantum number which occurs `count` times.
"""
QuantumNumbers(qn::AbstractQuantumNumber,count::Int=1)=QuantumNumbers('C',[qn],[0,count],qnsindptr)

"""
    QuantumNumbers(od::OrderedDict{<:AbstractQuantumNumber,Int})

Construct a `QuantumNumbers` from an ordered dict containing concrete quantum numbers and their counts.
"""
function QuantumNumbers(od::OrderedDict{<:AbstractQuantumNumber,Int})
    contents=Vector{od|>keytype}(undef,length(od))
    indptr=zeros(Int,length(od)+1)
    for (i,(qn,count)) in enumerate(od)
        @inbounds contents[i]=qn
        @inbounds indptr[i+1]=indptr[i]+count
    end
    QuantumNumbers('U',contents,indptr,qnsindptr)
end

"""
    QuantumNumbers(od::OrderedDict{<:AbstractQuantumNumber,UnitRange{Int}})

Construct a `QuantumNumbers` from an ordered dict containing concrete quantum numbers and their slices.
"""
function QuantumNumbers(od::OrderedDict{<:AbstractQuantumNumber,UnitRange{Int}})
    contents=Vector{od|>keytype}(undef,length(od))
    indptr=zeros(Int,length(od)+1)
    for (i,(qn,slice)) in enumerate(od)
        @inbounds contents[i]=qn
        @inbounds @assert indptr[i]+1==slice.start "QuantumNumbers error: slice not consistent."
        @inbounds indptr[i+1]=slice.stop
    end
    QuantumNumbers('U',contents,indptr,qnsindptr)
end

"""
    ==(qns1::QuantumNumbers,qns2::QuantumNumbers) -> Bool

Overloaded `==` operator. Two `QuantumNumbers`es are equal to each other if and only if both their `contents`es and `indptr`s are elementwise equal to each other.
!!! note
    It is not necessary for two `QuantumNumbers`es to have the same eltype nor the same form to be equal to each other.
"""
Base.:(==)(qns1::QuantumNumbers,qns2::QuantumNumbers)=all(qns1.contents.==qns2.contents) && all(qns1.indptr.==qns2.indptr)

"""
    show(io::IO,qns::QuantumNumbers)

Show a `QuantumNumbers`.
"""
Base.show(io::IO,qns::QuantumNumbers)=@printf io "QNS(%s)" join((@sprintf("%s=>%s",qn,slice) for (qn,slice) in pairs(qns,qnsindptr)),',')

"""
    string(qns::QuantumNumbers) -> String

Convert a `QuantumNumbers` to string.
"""
Base.string(qns::QuantumNumbers)=@sprintf "QNS(%s,%s)" qns|>length qns|>dimension

"""
    length(qns::QuantumNumbers) -> Int

Get the number of unduplicate qunatum numbers in the `QuantumNumbers`.
"""
Base.length(qns::QuantumNumbers)=length(qns.contents)

"""
    eltype(::Type{<:QuantumNumbers{QN}}) where QN
    eltype(qns::QuantumNumbers)

Get the type of the concrete `AbstractQuantumNumber` contained in a `QuantumNumbers`.
"""
Base.eltype(::Type{<:QuantumNumbers{QN}}) where QN=QN
Base.eltype(qns::QuantumNumbers)=qns|>typeof|>eltype

"""
    getindex(qns::QuantumNumbers,index::Int) -> AbstractQuantumNumber
    getindex(qns::QuantumNumbers,slice::UnitRange{Int}) -> QuantumNumbers
    getindex(qns::QuantumNumbers,indices::Vector{Int}) -> QuantumNumbers

Overloaded `[]` operator.
!!! note
    1. For a `QuantumNumbers`, all these `getindex` functions act on its `contents`, i.e. its compressed data, but not on its expansion, i.e. the uncompressed data. This definition is consistent with the [`length`](@ref) function.
    2. When the index is an integer, the result is a `AbstractQuantumNumber`, while when the index is a unit range or a vector of intgers, the result is a `QuantumNumbers`. The logic is quite reasonable because such behaviors are much alike to those of a vector container.
"""
Base.getindex(qns::QuantumNumbers,index::Int)=qns.contents[index]
function Base.getindex(qns::QuantumNumbers,slice::UnitRange{Int})
    contents=qns.contents[slice]
    indptr=qns.indptr[slice.start:slice.stop+1]
    indptr.=indptr.-qns.indptr[slice.start]
    return QuantumNumbers(qns.form,contents,indptr,qnsindptr)
end
function Base.getindex(qns::QuantumNumbers,indices::Vector{Int})
    contents=Vector{qns|>eltype}(undef,length(indices))
    indptr=zeros(Int,length(indices)+1)
    for (i,index) in enumerate(indices)
        contents[i]=qns.contents[index]
        indptr[i+1]=indptr[i]+qns.indptr[index+1]-qns.indptr[index]
    end
    return QuantumNumbers('U',contents,indptr,qnsindptr)
end

"""
    iterate(qns::QuantumNumbers,state::Int=1)
    iterate(rv::Iterators.Reverse{<:QuantumNumbers},state::Int=length(rv.itr,false))

Iterate or reversely iterate over the concrete `AbstractQuantumNumber`s contained in a `QuantumNumbers`.
"""
Base.iterate(qns::QuantumNumbers,state::Int=1)=state>length(qns) ? nothing : (@inbounds(qns.contents[state]),state+1)
Base.iterate(rv::Reverse{<:QuantumNumbers},state::Int=length(rv.itr))=state<1 ? nothing : (@inbounds(rv.itr.contents[state]),state-1)

"""
    keys(qns::QuantumNumbers) -> Vector{qns|>eltype}

Iterate over the concrete `AbstractQuantumNumber`s contained in a `QuantumNumbers`.
"""
Base.keys(qns::QuantumNumbers)=qns.contents

"""
    values(qns::QuantumNumbers,::typeof(qnsindptr))
    values(qns::QuantumNumbers,::typeof(qnscounts))

Iterate over the slices/counts of the `QuantumNumbers`.
"""
@views Base.values(qns::QuantumNumbers,::typeof(qnsindptr))=((start+1):stop for (start,stop) in zip(qns.indptr[1:end-1],qns.indptr[2:end]))
@views Base.values(qns::QuantumNumbers,::typeof(qnscounts))=(stop-start for (start,stop) in zip(qns.indptr[1:end-1],qns.indptr[2:end]))

"""
    pairs(qns::QuantumNumbers,choice::Union{typeof(qnsindptr),typeof(qnscounts)})

Iterate over the `AbstractQuantumNumber=>slice` or `AbstractQuantumNumber=>count` pairs.
"""
Base.pairs(qns::QuantumNumbers,choice::Union{typeof(qnsindptr),typeof(qnscounts)})=Base.Generator(=>,keys(qns),values(qns,choice))

"""
    dimension(qns::QuantumNumbers) -> Int

The dimension of the Hilbert space a `QuantumNumbers` represents.
"""
dimension(qns::QuantumNumbers)=@inbounds(qns.indptr[end])

"""
    sort(qns::QuantumNumbers) -> QuantumNumbers,Vector{Int}

Sort the quantum numbers of a `AbstractQuantumNumber`, return the sorted `AbstractQuantumNumber` and the permutation array that sorts the expansion of the original `QuantumNumbers`.
"""
function Base.sort(qns::QuantumNumbers)
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
    return QuantumNumbers('C',contents,indptr,qnsindptr),permutation
end

"""
    findall(qns::QuantumNumbers{QN},target::QN,choice::Union{typeof(qnscontents),typeof(qnsexpansion)}) where QN<:AbstractQuantumNumber -> Vector{Int}
    findall(qns::QuantumNumbers{QN},targets::NTuple{N,QN},::typeof(qnscontents)) where {N,QN<:AbstractQuantumNumber} -> Vector{Int}
    findall(qns::QuantumNumbers{QN},targets::NTuple{N,QN},::typeof(qnsexpansion)) where {N,QN<:AbstractQuantumNumber} -> Vector{Int}

Find all the indices of the target quantum numbers in the contents (`qnscontents` case) or the expansion (`qnsexpansion` case) of a `QuantumNumbers`.
"""
Base.findall(qns::QuantumNumbers{QN},target::QN,choice::Union{typeof(qnscontents),typeof(qnsexpansion)}) where QN<:AbstractQuantumNumber=findall(qns,(target,),choice)
function Base.findall(qns::QuantumNumbers{QN},targets::NTuple{N,QN},::typeof(qnscontents)) where {N,QN<:AbstractQuantumNumber}
    result=Int[]
    if qns.form=='C'
        for qn in targets
            range=searchsorted(qns.contents,qn)
            range.start<=range.stop && push!(result,range.start)
        end
    else
        for qn in targets
            append!(result,findall(isequal(qn),qns.contents))
        end
    end
    result
end
function Base.findall(qns::QuantumNumbers{QN},targets::NTuple{N,QN},::typeof(qnsexpansion)) where {N,QN<:AbstractQuantumNumber}
    result=Int[]
    for index in findall(qns,targets,qnscontents)
        @inbounds append!(result,(qns.indptr[index]+1):qns.indptr[index+1])
    end
    result
end

"""
    +(qn::AbstractQuantumNumber) -> AbstractQuantumNumber
    +(qn::QN,qns::QN...) where QN<:AbstractQuantumNumber -> QN
    +(qns::QuantumNumbers) -> QuantumNumbers
    +(qn::QN,qns::QuantumNumbers{QN}) where QN<:AbstractQuantumNumber -> QuantumNumbers{QN}
    +(qns::QuantumNumbers{QN},qn::QN) where QN<:AbstractQuantumNumber -> QuantumNumbers{QN}

Overloaded `+` operator for `AbstractQuantumNumber` and `QuantumNumbers`.
!!! note
    1. The addition between a `QuantumNumbers` and a `AbstractQuantumNumber` is just a global shift of the contents of the `QuantumNumbers` by the `AbstractQuantumNumber`, therefore, the result is a `QuantumNumbers`.
    2. `+` cannot be used between two `QuantumNumbers` because the result is ambiguous. Instead, use `⊕` for direct sum and `⊗` for direct product.
    3. To ensure type stability, two `AbstractQuantumNumber` can be added together if and only if they are of the same type.
    4. Similarly, a `AbstractQuantumNumber` and a `QuantumNumbers` can be added together if and only if the former's type is the same with the latter's eltype.
"""
Base.:+(qn::AbstractQuantumNumber)=qn
Base.:+(qn::QN,qns::QN...) where QN<:AbstractQuantumNumber=map(+,qn,qns...)
Base.:+(qns::QuantumNumbers)=qns
Base.:+(qn::QN,qns::QuantumNumbers{QN}) where QN<:AbstractQuantumNumber=qns+qn
Base.:+(qns::QuantumNumbers{QN},qn::QN) where QN<:AbstractQuantumNumber=QuantumNumbers(qns.form=='G' ? 'G' : 'U',[iqn+qn for iqn in qns],qns.indptr,qnsindptr)

"""
    -(qn::AbstractQuantumNumber) -> AbstractQuantumNumber
    -(qn1::QN,qn2::QN) where QN<:AbstractQuantumNumber -> QN
    -(qns::QuantumNumbers) -> QuantumNumbers
    -(qn::QN,qns::QuantumNumbers{QN}) where QN<:AbstractQuantumNumber -> QuantumNumbers{QN}
    -(qns::QuantumNumbers{QN},qn::QN) where QN<:AbstractQuantumNumber -> QuantumNumbers{QN}

Overloaded `-` operator for `AbstractQuantumNumber` and `QuantumNumbers`.
!!! note
    1. The subtraction between a `QuantumNumbers` and a `AbstractQuantumNumber` is just a global shift of the contents of the `QuantumNumbers` by the `AbstractQuantumNumber`, therefore, the result is a `QuantumNumbers`.
    2. `-` cannot be used between two `QuantumNumbers` because the result is ambiguous. Instead, use `⊕` with signs for direct sum and `⊗` with signs for direct product.
    3. To ensure type stability, a `AbstractQuantumNumber` can be subtracted by another `AbstractQuantumNumber` if and only if they are of the same type.
    4. Similarly, a `AbstractQuantumNumber` can be subtracted by a `QuantumNumbers` or vice versa if and only if the former's type is the same with the latter's eltype.
"""
Base.:-(qn::AbstractQuantumNumber)=map(-,qn)
Base.:-(qn1::QN,qn2::QN) where QN<:AbstractQuantumNumber=map(-,qn1,qn2)
Base.:-(qns::QuantumNumbers)=QuantumNumbers(qns.form=='G' ? 'G' : 'U',-qns.contents,qns.indptr,qnsindptr)
Base.:-(qn::QN,qns::QuantumNumbers{QN}) where QN<:AbstractQuantumNumber=QuantumNumbers(qns.form=='G' ? 'G' : 'U',[qn-iqn for iqn in qns],qns.indptr,qnsindptr)
Base.:-(qns::QuantumNumbers{QN},qn::QN) where QN<:AbstractQuantumNumber=QuantumNumbers(qns.form=='G' ? 'G' : 'U',[iqn-qn for iqn in qns],qns.indptr,qnsindptr)

"""
    *(qn::AbstractQuantumNumber,factor::Integer) -> AbstractQuantumNumber
    *(factor::Integer,qn::AbstractQuantumNumber) -> AbstractQuantumNumber
    *(qns::QuantumNumbers,factor::Integer) -> QuantumNumbers
    *(factor::Integer,qns::QuantumNumbers) -> QuantumNumbers

Overloaded `*` operator for the multiplication between an integer and a `AbstractQuantumNumber` or a `QuantumNumbers`.
"""
@generated function Base.:*(qn::AbstractQuantumNumber,factor::Integer)
    exprs=Expr[:(getfield(qn,$i)*factor) for i=1:(qn|>fieldnames|>length)]
    return :(typeof(qn)($(exprs...)))
end
Base.:*(factor::Integer,qn::AbstractQuantumNumber)=qn*factor
Base.:*(qns::QuantumNumbers,factor::Integer)=QuantumNumbers(qns.form=='G' ? 'G' : 'U',[qn*factor for qn in qns.contents],qns.indptr,qnsindptr)
Base.:*(factor::Integer,qns::QuantumNumbers)=qns*factor

"""
    ^(qn::AbstractQuantumNumber,factor::Integer) -> AbstractQuantumNumber
    ^(qns::QuantumNumbers,factor::Integer) -> QuantumNumbers

Overloaded `^` operator for `AbstractQuantumNumber` and `QuantumNumbers`. This operation translates into the direct product `⊗` of `factor` copies of `qn` or `qns`.
"""
Base.:^(qn::AbstractQuantumNumber,factor::Integer)=⊗(NTuple{factor,qn|>typeof}(qn for i=1:factor))
Base.:^(qns::QuantumNumbers,factor::Integer)=⊗(NTuple{factor,qns|>typeof}(qns for i=1:factor))

"""
    ⊕(qns::NTuple{N,<:AbstractQuantumNumber},signs::NTuple{N,Int}=ntuple(i->1,N)) where N -> QuantumNumbers
    ⊕(qnses::NTuple{N,QuantumNumbers{QN}},signs::NTuple{N,Int}=ntuple(i->1,N)) where {N,QN<:AbstractQuantumNumber} -> QuantumNumbers{QN}

Get the direct sum of some `AbstractQuantumNumber`s or `QuantumNumbers`es.
!!! note
    1. Physically, the direct sum of a couple of `AbstractQuantumNumber`s or `QuantumNumbers`es is defined by the direct sum of the bases of the Hilbert spaces they represent. Therefore, the input `AbstractQuantumNumber`s or `QuantumNumbers`es must be homogenous. Inhomogenous 'AbstractQuantumNumber's must be direct producted first to ensure homogenity before the direct sum.
    2. Apparently, the dimension of the result equals the summation of those of the inputs, which means, even for `AbstractQuantumNumber`s, the result will be naturally a `QuantumNumbers` because the dimension of the result is largeer than 1.
    3. Signs of `AbstractQuantumNumber`s or `QuantumNumbers`es can be provided when getting their direct sums.
"""
function ⊕(qns::NTuple{N,<:AbstractQuantumNumber},signs::NTuple{N,Int}=ntuple(i->1,N)) where N
    QuantumNumbers('G',[sign==1 ? qn : -qn for (sign,qn) in zip(signs,qns)],fill(1,qns|>length),qnscounts)
end
function ⊕(qnses::NTuple{N,QuantumNumbers{QN}},signs::NTuple{N,Int}=ntuple(i->1,N)) where {N,QN<:AbstractQuantumNumber}
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
    QuantumNumbers('G',contents,indptr,qnsindptr)
end

"""
    ⊗(::Type{QN},qn1::AbstractQuantumNumber,qn2::AbstractQuantumNumber) where QN<:AbstractQuantumNumber -> QN
    ⊗(qns::NTuple{N,<:AbstractQuantumNumber},signs::NTuple{N,Int}=ntuple(i->1,N)) where N -> AbstractQuantumNumber
    ⊗(qnses::NTuple{N,QuantumNumbers{QN}},signs::NTuple{N,Int}=ntuple(i->1,N)) where {N,QN<:AbstractQuantumNumber} -> QuantumNumbers{QN}

Get the direct product of some `AbstractQuantumNumber`s or `QuantumNumbers`es.
!!! note
    1. Physically, the direct product of a couple of `AbstractQuantumNumber`s or `QuantumNumbers`es are defined by the direct product of the bases of the Hilbert spaces they represent. Therefore, `QuantumNumbers` with differenct types or `QuantumNumbers`es with differenct eltypes are allowed to be direct producted in principle. However, for simplicity, we only implement a method which handle the situation of two `AbstractQuantumNumber`s with differenct types. The type of the result should be provided as the first parameter. Note that in this situation, the `fieldnames` and `periods` of the result type must be exactly equal to the flattened fieldnames and periods of the two input `AbstractQuantumNumber`s, which means, even the order of the input `AbstractQuantumNumber`s matters.
    2. Apparently, the dimension of the result equals the product of those of the inputs. Therefore, the direct product of `AbstractQuantumNumber`s is also a `AbstractQuantumNumber` since its dimension is still one.
    3. For other situations except the one mentioned in Note.1, the input `AbstractQuantumNumber`s or `QuantumNumbers`es must be homogenous. Meanwhile, signs can also be provided for these situations. Note that each quantum number in the contents of the result is obtained by a summation of the corresponding quanum numbers out of the inputs with the correct signs. This is a direct observation of the Abelian nature of our quantum numbers.
"""
function ⊗(::Type{QN},qn1::AbstractQuantumNumber,qn2::AbstractQuantumNumber) where QN<:AbstractQuantumNumber
    qnnames,qn1names,qn2names=QN|>fieldnames,qn1|>typeof|>fieldnames,qn2|>typeof|>fieldnames
    qnperiods,qn1periods,qn2periods=QN|>periods,qn1|>typeof|>periods,qn2|>typeof|>periods
    @assert qnnames==NTuple{QN|>length,Symbol}(flatten((qn1names,qn2names))) "⊗ error: fieldnames not match ($(qnnames),$(qn1names),$(qn2names))."
    @assert qnperiods==NTuple{QN|>length,Float64}(flatten((qn1periods,qn2periods))) "⊗ error: periods not match ($(qnperiods),$(qn1periods),$(qn2periods))."
    QN(qn1...,qn2...)
end
⊗(qns::NTuple{N,<:AbstractQuantumNumber},signs::NTuple{N,Int}=ntuple(i->1,N)) where N=sum(sign==1 ? qn : -qn for (sign,qn) in zip(signs,qns))
function ⊗(qnses::NTuple{N,QuantumNumbers{QN}},signs::NTuple{N,Int}=ntuple(i->1,N)) where {N,QN<:AbstractQuantumNumber}
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
    QuantumNumbers('G',contents,indptr,qnsindptr)
end

"""
    kron(::Type{QN},qn1::AbstractQuantumNumber,qn2::AbstractQuantumNumber) where QN<:AbstractQuantumNumber -> QN
    kron(qns::NTuple{N,<:AbstractQuantumNumber},signs::NTuple{N,Int}=ntuple(i->1,N)) where N -> AbstractQuantumNumber
    kron(qnses::NTuple{N,QuantumNumbers{QN}},signs::NTuple{N,Int}=ntuple(i->1,N)) where {N,QN<:AbstractQuantumNumber} -> QuantumNumbers{QN}

Kronecker product of some `AbstractQuantumNumber`s or `QuantumNumbers`es. This is defined to be equivalent to the direct product `⊗`.
"""
Base.kron(::Type{QN},qn1::AbstractQuantumNumber,qn2::AbstractQuantumNumber) where QN<:AbstractQuantumNumber=⊗(QN,qn1,qn2)
Base.kron(qns::NTuple{N,<:AbstractQuantumNumber},signs::NTuple{N,Int}=ntuple(i->1,N)) where N=⊗(qns,signs)
Base.kron(qnses::NTuple{N,QuantumNumbers{QN}},signs::NTuple{N,Int}=ntuple(i->1,N)) where {N,QN<:AbstractQuantumNumber}=⊗(qnses,signs)

"""
    ukron(qnses::NTuple{N,QuantumNumbers{QN}},signs::NTuple{N,Int}=ntuple(i->1,N)) where {N,QN<:AbstractQuantumNumber} -> QuantumNumbers{QN},Dict{QN,Dict{NTuple{N,QN},UnitRange{Int}}}

Unitary Kronecker product of several `QuantumNumbers`es. The product result as well as the records of the product will be returned.
!!! note
    1. All input `QuantumNumbers` must be 'U' formed or 'C' formed.
    2. Since duplicate quantum number are not allowed in 'U' formed and 'C' formed `QuantumNumbers`es, in general, there exists a merge process of duplicate quantum numbers in the product result. Therefore, records are needed to keep track of this process, which will be returned along with the product result. The records are stored in a `Dict{QN,Dict{NTuple{N,QN},UnitRange{Int}}}` typed dict, in which, for each unduplicate quantum number `qn` in the product result, there exist a record `Dict((qn₁,qn₂,...)=>start:stop,...)` telling what quantum numbers `(qn₁,qn₂,...)` a mereged duplicate `qn` comes from and what slice `start:stop` this merged duplicate corresponds.
"""
function ukron(qnses::NTuple{N,QuantumNumbers{QN}},signs::NTuple{N,Int}=ntuple(i->1,N)) where {N,QN<:AbstractQuantumNumber}
    @assert all(qns.form=='U' || qns.form=='C' for qns in qnses) "⊗ error: all input qnses should be 'U' formed or 'C' formed."
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
    QuantumNumbers('C',contents,indptr,qnsindptr),records
end

"""
    expand(qns::QuantumNumbers,::typeof(qnscontents)) -> Vector{qns|>eltype}
    expand(qns::QuantumNumbers,::typeof(qnsindices)) -> Vector{Int}

Expand the contents (`qnscontents` case) or indices (`qnsindices` case) of a `QuantumNumbers` to the uncompressed form.
"""
function expand(qns::QuantumNumbers,::typeof(qnscontents))
    result=Vector{qns|>eltype}(undef,dimension(qns))
    for i=1:length(qns)
        for j=@inbounds((qns.indptr[i]+1):qns.indptr[i+1])
            @inbounds result[j]=qns.contents[i]
        end
    end
    return result
end
function expand(qns::QuantumNumbers,::typeof(qnsindices))
    result=Vector{Int}(undef,dimension(qns))
    for i=1:length(qns)
        for j=@inbounds((qns.indptr[i]+1):qns.indptr[i+1])
            @inbounds result[j]=i
        end
    end
    return result
end

"""
    decompose(qnses::NTuple{N,QuantumNumbers{QN}},target::QN,signs::NTuple{N,Int},::typeof(qnsbruteforce);nmax::Int=20) where {N,QN<:AbstractQuantumNumber} -> Vector{NTuple{N,Int}}
    decompose(qnses::NTuple{N,QuantumNumbers{QN}},target::QN,signs::NTuple{N,Int},::typeof(qnsmontecarlo);nmax::Int=20) where {N,QN<:AbstractQuantumNumber} -> Vector{NTuple{N,Int}}

Find a couple of decompositions of `target` with respect to `qnses`.
!!! note
    A tuple of integers `(i₁,i₂,...)` is called a decomposition of a given `target` with respect to the given `qnses` if and only if they satisfy the "decomposition rule":
    ```math
    \\sum_\\text{j} \\text{signs}[\\text{j}]\\times\\text{qnses}[\\text{j}][\\text{i}_{\\text{j}}]==\\text{target}
    ```
    This equation is in fact a kind of a set of restricted [linear Diophantine equations](https://en.wikipedia.org/wiki/Diophantine_equation#Linear_Diophantine_equations). Indeed, our quantum numbers are always discrete Abelian ones and all instances of a concrete `AbstractQuantumNumber` forms a [module](https://en.wikipedia.org/wiki/Module_(mathematics)) over the [ring](https://en.wikipedia.org/wiki/Ring_(mathematics)) of integers. Therefore, each quantum number can be represented as a integral multiple of the unit element of the Abelian module, which results in the final reduction of the above equation to a set of linear Diophantine equations. Then finding a decomposition is equivalent to find a solution of the reduced linear Diophantine equations, with the restriction that the quantum numbers constructed from the solution should be in the corresponding `qnses`. Here we provide two methods to find such decompositions, one is by brute force (`qnsbruteforce` case), and the other is by Monte Carlo simultatioins (`qnsmontecarlo` case).
"""
function decompose(qnses::NTuple{N,QuantumNumbers{QN}},target::QN,signs::NTuple{N,Int},::typeof(qnsbruteforce);nmax::Int=20) where {N,QN<:AbstractQuantumNumber}
    result=Set{NTuple{N,Int}}()
    cache=Vector{Int}(undef,N)
    dimensions=NTuple{N,Int}(dimension(qns) for qns in reverse(qnses))
    indices=findall(⊗(qnses,signs),(target,),qnsexpansion)
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
function decompose(qnses::NTuple{N,QuantumNumbers{QN}},target::QN,signs::NTuple{N,Int},::typeof(qnsmontecarlo);nmax::Int=20) where {N,QN<:AbstractQuantumNumber}
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
    subset(qns::QuantumNumbers{QN},target::QN) where QN<:AbstractQuantumNumber -> QuantumNumbers{QN}
    subset(qns::QuantumNumbers{QN},targets::NTuple{N,QN}) where {N,QN<:AbstractQuantumNumber} -> QuantumNumbers{QN}

Find a subset of a `QuantumNumbers` by picking out the quantum numbers in targets.
"""
subset(qns::QuantumNumbers{QN},target::QN) where QN<:AbstractQuantumNumber=qns[findall(qns,target,qnscontents)]
subset(qns::QuantumNumbers{QN},targets::NTuple{N,QN}) where {N,QN<:AbstractQuantumNumber}=qns[findall(qns,targets,qnscontents)]

"""
    reorder(qns::QuantumNumbers,permutation::Vector{Int},::typeof(qnscontents)) -> QuantumNumbers
    reorder(qns::QuantumNumbers,permutation::Vector{Int},::typeof(qnsexpansion)) -> QuantumNumbers

Reorder the quantum numbers contained in a `QuantumNumbers` with a permutation and return the new one. For `qnscontents` case, the permutation is for the contents of the original `QuantumNumbers` while for `qnsexpansion` case, the permutation is for the expansion of the original `QuantumNumbers`.
"""
reorder(qns::QuantumNumbers,permutation::Vector{Int},::typeof(qnscontents))=qns[permutation]
function reorder(qns::QuantumNumbers,permutation::Vector{Int},::typeof(qnsexpansion))
    contents=Vector{qns|>eltype}(undef,length(permutation))
    indptr=zeros(Int,length(permutation)+1)
    expansion=expand(qns,qnsindices)
    for (i,p) in enumerate(permutation)
        contents[i]=qns.contents[expansion[p]]
        indptr[i+1]=indptr[i]+1
    end
    return QuantumNumbers('G',contents,indptr,qnsindptr)
end

"""
    toordereddict(qns::QuantumNumbers,::typeof(qnsindptr)) -> OrderedDict{qns|>eltype,UnitRange{Int}}
    toordereddict(qns::QuantumNumbers,::typeof(qnscounts)) -> OrderedDict{qns|>eltype,Int}

Convert a `QuantumNumbers` to an ordered dict.
"""
function toordereddict(qns::QuantumNumbers,::typeof(qnsindptr))
    @assert qns.form != 'G' "toordereddict error: input `QuantumNumbers` cannot be `G` formed."
    result=OrderedDict{qns|>eltype,UnitRange{Int}}()
    for i=1:length(qns)
        @inbounds result[qns.contents[i]]=qns.indptr[i]+1:qns.indptr[i+1]
    end
    return result
end
function toordereddict(qns::QuantumNumbers,::typeof(qnscounts))
    @assert qns.form != 'G' "toordereddict error: input `QuantumNumbers` cannot be `G` formed."
    result=OrderedDict{qns|>eltype,Int}()
    for i=1:length(qns)
        @inbounds result[qns.contents[i]]=qns.indptr[i+1]-qns.indptr[i]
    end
    return result
end

"""
    SQNS(S::Real)

Construct the `QuantumNumbers` of the Hilbert space of a signle spin `S`.
"""
SQNS(S::Real)=QuantumNumbers('C',[SQN(sz) for sz=-S:S],collect(0:Int(2*S+1)),qnsindptr)

"""
    PQNS(N::Real)

Construct the `QuantumNumbers` of the Hilbert space of a single-particle state with at most `N` identical particles.
"""
PQNS(N::Real)=QuantumNumbers('C',[PQN(np) for np=0:N],collect(0:Int(N)+1),qnsindptr)

"""
    SzPQNS(Sz::Real)

Construct the `QuantumNumbers` of the Hilbert space of a single-paritcle state with at most one particle whose spin-z component is `Sz`.
"""
SzPQNS(Sz::Real)=QuantumNumbers('C',[SPQN(0.0,0.0),SPQN(1.0,Sz)],[0,1,2],qnsindptr)

"""
    SPQNS(S::Real)

Construct the `QuantumNumbers` of the Hilbert space of a single site with internal degrees of freedom that can be ascribed to a spin `S`.
"""
function SPQNS(S::Real)
    sqns=[SQN(sz) for sz=-S:S]
    contents=[SPQN(0.0,0.0)]
    for n=1:length(sqns)
        pn=PQN(n*1.0)
        for spins in combinations(sqns,n)
            push!(contents,⊗(SPQN,pn,sum(spins)))
        end
    end
    return sort(QuantumNumbers('G',contents,collect(0:length(contents)),qnsindptr))[1]
end

"""
    Z2QNS()

Construct the `QuantumNumbers` of a ``Z_2`` Hilbert space.
"""
Z2QNS()=QuantumNumbers('C',[Z2QN(0.0),Z2QN(1.0)],[0,1,2],qnsindptr)

end #module
