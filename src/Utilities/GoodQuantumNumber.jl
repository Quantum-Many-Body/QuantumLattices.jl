module GoodQuantumNumber

import Printf: @printf,@sprintf
import Base.Enums: @enum
import ..NamedVector: AbstractNamedVector

export QuantumNumber
export @quantumnumber,periods,SQN,PQN,SPQN,Z2QN
export QNSProtocol,qnscounts,qnsindptr
export QuantumNumbers
export regularize!,regularization,⊕,⊗

"Abstract type for all concrete quantum numbers for a single basis."
abstract type QuantumNumber{A<:AbstractVector{Float64}} <: AbstractNamedVector{Float64,A} end

"""
    @quantumnumber typename fieldnames fieldperiods

Construct a concrete `QuantumNumber` with the type name being `typename`, fieldnames specified by `fieldnames` and periods specified by `fieldperiods`.
"""
macro quantumnumber(typename,fieldnames,fieldperiods)
    typename=Symbol(typename)
    fieldnames=tuple(eval(fieldnames)...)
    fieldperiods=tuple(eval(fieldperiods)...)
    @assert length(fieldnames)==length(fieldperiods) "quantumnumber error: number of fieldnames($(length(fieldnames))) and fieldperiods($(length(fieldperiods))) not equal."
    @assert all(isa(name,Symbol) for name in fieldnames) "quantumnumber error: all field names should be Symbol."
    @assert all(fieldperiods.>0) "quantumnumber error: all field periods should be greater than 0."
    return quote
        global periods
        struct $(esc(typename)){A<:AbstractVector{Float64}} <: QuantumNumber{A}
            values::A
            $(esc(typename))(values::AbstractVector{Float64},regularize::Bool=true)=(regularize && regularize!($(esc(typename)),values);new{values|>typeof}(values))
            $(esc(typename))(values::Real...)=$(esc(typename))(collect(Float64,values))
        end
        $(esc(typename)){A}(values::AbstractVector{Float64},regularize::Bool=true) where A=$(esc(typename))(values,regularize)
        Base.fieldnames(::Type{<:$(esc(typename))},private=false)=private ? tuple(($fieldnames)...,:values) : $fieldnames
        periods(::Type{<:$(esc(typename))})=$fieldperiods
    end
end

"""
    SQN(Sz::Real)

The concrete `QuantumNumber` of a quantum system with spin z-component `Sz` conserved.
"""
SQN
@quantumnumber "SQN" (:Sz,) (Inf,)

"""
    PQN(N::Real)

The concrete `QuantumNumber` of a quantum system with particle number `N` conserved.
"""
PQN
@quantumnumber "PQN" (:N,) (Inf,)

"""
    SPQN(N::Real,Sz::Real)

The concrete `QuantumNumber` of a quantum system with both particle number `N` and spin z-component `Sz` conserved.
"""
SPQN
@quantumnumber "SPQN" (:N,:Sz) (Inf,Inf)

"""
    Z2QN(N::Real)

The concrete `QuantumNumber` of a quantum system with a Z₂-like conserved quantity.
"""
Z2QN
@quantumnumber "Z2QN" (:N,) (2,)

"""
Protocol used for the initilization of `QuantumNumbers`.
* `qnscounts`: initilization by counts
* `qnsindptr`: initilization by indptr
"""
@enum QNSProtocol qnscounts=0 qnsindptr=1

"""
    QuantumNumbers(form::Char,::Type{QN},contents::Array{Float64,2},info::Array{Int,1},protocol::QNSProtocol=qnscounts,regularize::Bool=true) where QN<:QuantumNumber
    QuantumNumbers(form::Char,contents::Array{QN,1},info::Array{Int,1},protocol::QNSProtocol=qnscounts) where QN<:QuantumNumber
    QuantumNumbers(qn::QuantumNumber,count::Int=1)

The whole quantum numbers of the total bases of a Hilbert space.
"""
struct QuantumNumbers{QN<:QuantumNumber}
    form::Char
    contents::Array{Float64,2}
    indptr::Array{Int,1}
    function QuantumNumbers(form::Char,::Type{QN},contents::Array{Float64,2},info::Array{Int,1},protocol::QNSProtocol=qnscounts,regularize::Bool=true) where QN<:QuantumNumber
        @assert form|>uppercase ∈ Set(['G','U','C']) "QuantumNumbers error: 'form'($form) is not 'G','U' or 'C'."
        @assert protocol ∈ Set([qnscounts,qnsindptr]) "QuantumNumbers error: protocol($protocol) not supported."
        regularize!(QN,contents)
        csize,isize=size(contents,2),size(info,1)
        if protocol==qnscounts
            @assert csize==isize "QuantumNumbers error: dismatch shapes of contents and info ($(csize)!=$isize)."
            indptr=[0,cumsum(info)...]
        else
            @assert csize+1==isize "QuantumNumbers error: dismatch shapes of contents and info ($csize+1!=$isize)."
            indptr=info
        end
        new{QN}(form|>uppercase,contents,indptr)
    end
    function QuantumNumbers(form::Char,contents::Array{QN,1},info::Array{Int,1},protocol::QNSProtocol=qnscounts) where QN<:QuantumNumber
        QuantumNumbers(form,QN,hcat([[qn.values...] for qn in contents]...),info,protocol,false)
    end
    QuantumNumbers(qn::QuantumNumber,count::Int=1)=QuantumNumbers('C',[qn],[count],qnscounts)
end

"""
    show(io::IO,qns::QuantumNumbers)

Show a `QuantumNumbers`.
"""
Base.show(io::IO,qns::QuantumNumbers)=@printf io "QNS(%s)" join((@sprintf("%s=>%s",qn,slice) for (qn,slice) in pairs(qns)),',')

"""
    string(qns::QuantumNumbers)

Convert a `QuantumNumbers` to string.
"""
Base.string(qns::QuantumNumbers)=@sprintf "QNS(%s,%s)" length(qns,false) length(qns,true)

"""
    length(qns::QuantumNumbers,duplicate::Bool=true)

Get the number of qunatum numbers in the `QuantumNumbers`.
* `duplicate==true`: the duplicate quantum numbers are counted duplicately. Then the result equals the dimension of the `QuantumNumbers`.
* `duplicate==false`: only unduplicate quantum numbers are counted. Then the result equals the number of columns of the `QuantumNumbers`' `contents`.
"""
Base.length(qns::QuantumNumbers,duplicate::Bool=true)=duplicate ? qns.indptr[end] : size(qns.contents,2)

"""
    eltype(::Type{<:QuantumNumbers{QN}}) where QN

Get the type of the concrete `QuantumNumber` contained in `QuantumNumbers`.
"""
Base.eltype(::Type{<:QuantumNumbers{QN}}) where QN=QN

"""
    getindex(qns::QuantumNumbers,index::Int)

Overloaded `[]` operator.
"""
@views Base.getindex(qns::QuantumNumbers,index::Int)=(qns|>typeof|>eltype)(qns.contents[:,index],false)

"""
    iterate(qns::QuantumNumbers,state::Int=1)

Iterate over the concrete `QuantumNumber`s the `QuantumNumbers` contains.
"""
Base.iterate(qns::QuantumNumbers,state::Int=1)=state>length(qns,false) ? nothing : (qns[state],state+1)

"""
    keys(qns::QuantumNumbers)

Iterate over the concrete `QuantumNumber`s the `QuantumNumbers` contains.
"""
Base.keys(qns::QuantumNumbers)=(qns[i] for i in 1:length(qns,false))

"""
    values(qns::QuantumNumbers,protocol::QNSProtocol=qnsindptr)

Iterate over the slices/counts of the `QuantumNumbers`.
"""
@views Base.values(qns::QuantumNumbers,protocol::QNSProtocol=qnsindptr)=(protocol==qnsindptr ? (s+1:e) : e-s for (s,e) in zip(qns.indptr[1:end-1],qns.indptr[2:end]))

"""
    pairs(qns::QuantumNumbers,protocol::QNSProtocol=qnsindptr)

Iterate over the `QuantumNumber=>slice` or `QuantumNumber=>count` pairs.
"""
Base.pairs(qns::QuantumNumbers,protocol::QNSProtocol=qnsindptr)=Base.Generator(=>,keys(qns),values(qns,protocol))

"""
    regularize!(::Type{QN},array::AbstractVector{Float64}) where QN<:QuantumNumber
    regularize!(::Type{QN},array::AbstractArray{Float64,2}) where QN<:QuantumNumber

Regularize an array representation of a `QuantumNumber`/`QuantumNumbers` by the concrete `QuantumNumber`'s periods.
"""
regularize!

function regularize!(::Type{QN},array::AbstractVector{Float64}) where QN<:QuantumNumber
    @assert size(array,1)==QN|>periods|>length
    for (i,period) in enumerate(QN|>periods)
        period!=Inf && (array[i]=array[i]%period+(array[i]%period<0 ? period : 0))
    end
end

function regularize!(::Type{QN},array::AbstractArray{Float64,2}) where QN<:QuantumNumber
    @assert size(array,1)==QN|>periods|>length
    for (i,period) in enumerate(QN|>periods)
        if period!=Inf
            for j in 1:size(array,2)
                array[i,j]=array[i,j]%period+(array[i,j]%period<0 ? period : 0)
            end
        end
    end
end

"""
    regularization(::Type{QN},array::AbstractArray{Float64,N}) where {QN<:QuantumNumber,N}

Get the regularized array of the array representation of a `QuantumNumber`/`QuantumNumbers` by the periods.
"""
function regularization(::Type{QN},array::AbstractArray{Float64,N}) where {QN<:QuantumNumber,N}
    @assert N==1||N==2 "regularization error: the dimension($(N)) of the input array must be one or two."
    result=copy(array)
    regularize!(QN,result)
    return result
end

"""
    +(qn::QuantumNumber)
    +(qn1::QN,qn2::QN) where QN<:QuantumNumber
    +(qn1::QN,qn2::QN,qns::QN...) where QN<:QuantumNumber
    +(qns::QuantumNumbers)
    +(qn::QN,qns::QuantumNumbers{QN}) where QN
    +(qns::QuantumNumbers{QN},qn::QN) where QN

Overloaded `+` operator for `QuantumNumber` and `QuantumNumbers`.
"""
Base.:+

Base.:+(qn::QuantumNumber)=qn
Base.:+(qn1::QN,qn2::QN) where QN<:QuantumNumber=QN(getfield(qn1,:values)+getfield(qn2,:values))
Base.:+(qn1::QN,qn2::QN,qns::QN...) where QN<:QuantumNumber=QN([sum(vs) for vs in zip(qn1,qn2,qns...)])
Base.:+(qns::QuantumNumbers)=qns
Base.:+(qn::QN,qns::QuantumNumbers{QN}) where QN=qns+qn
Base.:+(qns::QuantumNumbers{QN},qn::QN) where QN=nothing

"""
    -(qn::QuantumNumber)
    -(qn1::QN,qn2::QN) where QN<:QuantumNumber
    -(qns::QuantumNumbers)
    -(qn::QN,qns::QuantumNumbers{QN}) where QN
    -(qns::QuantumNumbers{QN},qn::QN) where QN

Overloaded `-` operator for `QuantumNumber` and `QuantumNumbers`.
"""
Base.:-

Base.:-(qn::QuantumNumber)=typeof(qn)(-getfield(qn,:values))
Base.:-(qn1::QN,qn2::QN) where QN<:QuantumNumber=QN(getfield(qn1,:values)-getfield(qn2,:values))
Base.:-(qns::QuantumNumbers)=QuantumNumbers(qns.form,qns|>typeof|>eltype,-qns.contents,qns.indptr,qnsindptr)
Base.:-(qn::QN,qns::QuantumNumbers{QN}) where QN=nothing
Base.:-(qns::QuantumNumbers{QN},qn::QN) where QN=nothing

"""
    *(qn::QuantumNumber,factor::Integer)
    *(factor::Integer,qn::QuantumNumber)
    *(qns::QuantumNumbers,factor::Integer)
    *(factor::Integer,qns::QuantumNumbers)

Overloaded `*` operator for `QuantumNumber` and `QuantumNumbers`..
"""
Base.:*

Base.:*(qn::QuantumNumber,factor::Integer)=factor*qn
Base.:*(factor::Integer,qn::QuantumNumber)=typeof(qn)(getfield(qn,:values)*factor)
Base.:*(qns::QuantumNumbers,factor::Integer)=factor*qns
Base.:*(factor::Integer,qns::QuantumNumbers)=QuantumNumbers(qns.form,qns|>typeof|>eltype,qns.contents*factor,qns.indptr,qnsindptr)

"""
    ^(qn::QuantumNumber,factor::Integer)
    ^(qns::QuantumNumbers,factor::Integer)

Overloaded `^` operator for `QuantumNumber` and `QuantumNumbers`.
"""
Base.:^

Base.:^(qn::QuantumNumber,factor::Integer)=QuantumNumbers(qn,factor)
Base.:^(qns::QuantumNumbers,factor::Integer)=⊗(fill(qns,factor))

"""
    ⊕(qns::QN...;signs::Union{Array{Int,1},Nothing}=nothing) where QN<:QuantumNumber
    ⊕(qnses::QuantumNumbers{QN}...;signs::Union{Array{Int,1},Nothing}=nothing) where QN

Get the direct sum of some `QuantumNumber`s or `QuantumNumbers`s.
"""
⊕

⊕(qns::QN...;signs::Union{Array{Int,1},Nothing}=nothing) where QN<:QuantumNumber=nothing
⊕(qnses::QuantumNumbers{QN}...;signs::Union{Array{Int,1},Nothing}=nothing) where QN=nothing

"""
    ⊗(::Type{QN},qn1::QuantumNumber,qn2::QuantumNumber) where QN<:QuantumNumber
    ⊗(qnses::QuantumNumbers{QN}...;signs::Union{Array{Int,1},Nothing}=nothing) where QN

Get the direct product of some `QuantumNumber`s or `QuantumNumbers`s.
"""
⊗

function ⊗(::Type{QN},qn1::QuantumNumber,qn2::QuantumNumber) where QN<:QuantumNumber
    qnnames,qn1names,qn2names=QN|>fieldnames,qn1|>typeof|>fieldnames,qn2|>typeof|>fieldnames
    qnperiods,qn1periods,qn2periods=QN|>periods,qn1|>typeof|>periods,qn2|>typeof|>periods
    @assert qnnames==tuple(qn1names...,qn2names...) "⊗ error: fieldnames not match ($(qnnames),$(qn1names),$(qn2names))."
    @assert qnperiods==tuple(qn1periods...,qn2periods...) "⊗ error: periods not match ($(qnperiods),$(qn1periods),$(qn2periods))."
    QN([getfield(qn1,:values);getfield(qn2,:values)])
end

⊗(qnses::QuantumNumbers{QN}...;signs::Union{Array{Int,1},Nothing}=nothing) where QN=nothing

end #module
