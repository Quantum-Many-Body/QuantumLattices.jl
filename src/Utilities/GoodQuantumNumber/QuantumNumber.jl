import Printf: @printf

export QuantumNumber,regularization,⊕

"Abstract type for all concrete quantum numbers for a single basis."
abstract type QuantumNumber end

"""
    qn.key

Overloaded `.` operator.
"""
function Base.getproperty(qn::QuantumNumber,key::Symbol)
    if key==:values
        result=getfield(qn,key)
    elseif key ∈ typeof(qn)|>fieldnames
        result=getfield(qn,:values)[findfirst(isequal(key),typeof(qn)|>fieldnames)]
    else
        error("$(string(typeof(qn))) getproperty error: :$key not available.")
    end
    result
end

"""
    ==(qn1::T,qn2::T) where T<:QuantumNumber

Overloaded `==` operator.
"""
Base.:(==)(qn1::T,qn2::T) where T<:QuantumNumber=getfield(qn1,:values)==getfield(qn2,:values)

"""
    show(io::IO,qn::QuantumNumber)

Show a concrete `QuantumNumber`.
"""
Base.show(io::IO,qn::QuantumNumber)=@printf io "%s(%s)" typeof(qn) join(getfield(qn,:values),',')

"""
    hash(qn::QuantumNumber,h::UInt)

Hash a concrete `QuantumNumber`.
"""
Base.hash(qn::QuantumNumber,h::UInt)=hash(getfield(qn,:values),h)

"""
    length(qn::QuantumNumber)

Get the length of a concrete `QuantumNumber`.
"""
Base.length(qn::QuantumNumber)=length(getfield(qn,:values))

"""
    eltype(::Type{<:QuantumNumber})

Get the data type of a concrete `QuantumNumber`. This is defined to be `Float64` by default.
"""
Base.eltype(::Type{<:QuantumNumber})=Float64

"""
    zero(::Type{QN}) where QN<:QuantumNumber

Get a concrete `QuantumNumber` with all values being zero.
"""
Base.zero(::Type{QN}) where QN<:QuantumNumber=QN((zero(Float64) for i in 1:(QN|>periods|>length))...)

"""
    ⊕(::Type{QN},qn1::QuantumNumber,qn2::QuantumNumber) where {QN<:QuantumNumber}

Get the directsum of two concrete `QuantumNumber`s. The returned type is specified by `QN`.
"""
function ⊕(::Type{QN},qn1::QuantumNumber,qn2::QuantumNumber) where {QN<:QuantumNumber}
    qnnames,qn1names,qn2names=QN|>fieldnames,qn1|>typeof|>fieldnames,qn2|>typeof|>fieldnames
    qnperiods,qn1periods,qn2periods=QN|>periods,qn1|>typeof|>periods,qn2|>typeof|>periods
    @assert qnnames==tuple(qn1names...,qn2names...) "⊕ error: fieldnames not match ($(qnnames),$(qn1names),$(qn2names))."
    @assert qnperiods==tuple(qn1periods...,qn2periods...) "⊕ error: periods not match ($(qnperiods),$(qn1periods),$(qn2periods))."
    QN(tuple(getfield(qn1,:values)...,getfield(qn2,:values)...))
end

"""
    regularization(::Type{QN},array::Array{Float64,1}) where QN<:QuantumNumber

Regularize an array representation of a concrete `QuantumNumber` by its periods.
"""
function regularization(::Type{QN},array::Array{Float64,1}) where QN<:QuantumNumber
    @assert size(array,1)==QN|>periods|>length
    for (i,period) in enumerate(QN|>periods)
        period!==nothing && (array[i]÷=period)
    end
    array
end
