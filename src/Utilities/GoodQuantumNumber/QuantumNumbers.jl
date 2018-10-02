export qnscounts,qnsindptr,QuantumNumbers,regularization

"Protocol used for the initilization of `QuantumNumbers`."
const qnscounts=0

"Protocol used for the initilization of `QuantumNumbers`"
const qnsindptr=1

"""
    QuantumNumbers(form::Char,::Type{QN},contents::Array{Float64,2},info::Array{Int64,1},protocol::Int64=qnscounts) where QN<:QuantumNumber
    QuantumNumbers(form::Char,contents::Array{QN,2},info::Array{Int64,1},protocol::Int64=qnscounts) where QN<:QuantumNumber

The whole quantum numbers of the total bases of a Hilbert space.
"""
struct QuantumNumbers{QN<:QuantumNumber}
    form::Char
    contents::Array{Float64,2}
    indptr::Array{Int64,1}
    function QuantumNumbers(form::Char,::Type{QN},contents::Array{Float64,2},info::Array{Int64,1},protocol::Int64=qnscounts) where QN<:QuantumNumber
        @assert protocol ∈ Set([qnscounts,qnsindptr]) "QuantumNumbers error: protocol($protocol) not supported."
        @assert form|>uppercase ∈ Set(['G','U','C']) "QuantumNumbers error: 'form'($form) is not 'G','U' or 'C'."
        if protocol==qnscounts
            @assert size(contents,2)==size(info,1) "QuantumNumbers error: dismatch shapes of contents and info ($(size(contents,2))!=$(size(info,1)))."
            indptr=[0,cumsum(info)...]
        else
            @assert size(contents,2)+1==size(info,1) "QuantumNumbers error: dismatch shapes of contents and info ($(size(contents,2)+1)!=$(size(info,1)))."
            indptr=info
        end
        new{QN}(form|>uppercase,contents,indptr)
    end
end

QuantumNumbers(form::Char,contents::Array{QN,2},info::Array{Int64,1},protocol::Int64=qnscounts) where QN<:QuantumNumber=QuantumNumbers(form,QN,hcat([[qn.values...] for qn in contents]...),info,protocol)

"""
    eltype(::Type{<:QuantumNumbers{QN}}) where QN

Get the type of the concrete `QuantumNumber` contained in `QuantumNumbers`.
"""
eltype(::Type{<:QuantumNumbers{QN}}) where QN=QN

"""
    regularization(::Type{QuantumNumbers{QN}},array::Array{Float64,2}) where QN

Regularize an array representation of a `QuantumNumbers` by the periods of the concrete `QuantumNumber` it contains.
"""
function regularization(::Type{QuantumNumbers{QN}},array::Array{Float64,2}) where QN
    @assert size(array,1)==QN|>periods|>length
    for (i,period) in enumerate(QN|>periods)
        period!==nothing && (array[i,:]÷=array[i,:])
    end
    array
end
