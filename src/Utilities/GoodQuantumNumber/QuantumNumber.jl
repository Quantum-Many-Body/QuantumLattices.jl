"""
    QuantumNumber{T<:Qndtypes}

Abstract type for all concrete quantum numbers for a single basis.
"""
abstract type QuantumNumber{T<:Qndtypes} end

==(qn1::QuantumNumber{T},qn2::QuantumNumber{T}) where T=getfield(qn1,:values)==getfield(qn2,:values)
show(io::IO,qn::QuantumNumber)=@printf io "%s(%s)" typeof(qn) join(getfield(qn,:values),',')
hash(qn::QuantumNumber{T},h::UInt) where T=hash(getfield(qn,:values),h)

"The data type of the quantum number."
qndtype(::Type{<:QuantumNumber{T}}) where T=T


"""
    SPQN(n::T1,sz::T2) where {T1<:Qndtypes,T2<:Qndtypes}

The concrete `QuantumNumber` of a quantum system with both particle number and spin z-component conserved.
* function `N`: the particle number
* function `Sz`: the spin z-component
"""
struct SPQN <: QuantumNumber{Float64}
    values::Tuple{Float64,Float64}
    SPQN(n::T1,sz::T2) where {T1<:Qndtypes,T2<:Qndtypes}=new((Float64(n),Float64(sz)))
end

names(::Type{SPQN})=("N","Sz")
periods(::Type{SPQN})=(nothing,nothing)
N(qn::SPQN)=qn.values[1]
Sz(qn::SPQN)=qn.values[2]
