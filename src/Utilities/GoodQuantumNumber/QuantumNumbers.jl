const qnscounts=0
const qnsindptr=1

"""
    QuantumNumbers(form::Char,qntype::Type{QN},contents::Array{T,2},counts::Array{Int64,1},::Val{qnscounts}) where {T<:Qndtypes,QN<:QuantumNumber{T}})
    QuantumNumbers(form::Char,qntype::Type{QN},contents::Array{T,2},counts::Array{Int64,1},::Val{qnscounts}) where {T<:Qndtypes,QN<:QuantumNumber{T}}

The whole quantum numbers of the total base of a Hilbert space.
"""
struct QuantumNumbers{QN<:QuantumNumber,T<:Qndtypes}
    form::Char
    contents::Array{T,2}
    indptr::Array{Int64,1}
    function QuantumNumbers(form::Char,qntype::Type{QN},contents::Array{T,2},indptr::Array{Int64,1},::Val{qnsindptr}) where {T<:Qndtypes,QN<:QuantumNumber{T}}
        @assert form|>uppercase ∈ Set(['G','U','C']) "QuantumNumbers construction error: 'form'($form) is not 'G','U' or 'C'."
        @assert size(contents,2)+1==size(indptr,1) "QuantumNumbers construction error: size(contents,2)+1!=size((indptr,1) ($(size(contents,2)+1)!=$(size(indptr,1))."
        new{QN,T}(form|>uppercase,contents,indptr)
    end
end

function QuantumNumbers(form::Char,qntype::Type{QN},contents::Array{T,2},counts::Array{Int64,1},::Val{qnscounts}) where {T<:Qndtypes,QN<:QuantumNumber{T}}
    QuantumNumbers(form,qntype,contents,[0,cumsum(counts)...],Val(qnsindptr))
end

qndtype(::Type{<:QuantumNumbers{QN,T}}) where {QN,T}=T
qntype(::Type{<:QuantumNumbers{QN,T}}) where {QN,T}=QN
