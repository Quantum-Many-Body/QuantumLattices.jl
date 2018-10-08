export SPQN,periods

"""
    SPQN(n::T1,sz::T2) where {T1<:Real,T2<:Real}

The concrete `QuantumNumber` of a quantum system with both particle number `n` and spin z-component `sz` conserved.
"""
struct SPQN <: QuantumNumber
    values::NTuple{2,Float64}
    SPQN(n::T1,sz::T2) where {T1<:Real,T2<:Real}=new((Float64(n),Float64(sz)))
end

Base.fieldnames(::Type{SPQN},private=false)=private ? (:N,:Sz,:values) : (:N,:Sz)
periods(::Type{SPQN})=(nothing,nothing)
