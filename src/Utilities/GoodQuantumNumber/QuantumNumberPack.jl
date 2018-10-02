export SPQN,periods,N,Sz

"""
    SPQN(n::T1,sz::T2) where {T1<:Real,T2<:Real}

The concrete `QuantumNumber` of a quantum system with both particle number `n` and spin z-component `sz` conserved.
"""
struct SPQN <: QuantumNumber
    values::NTuple{2,Float64}
    SPQN(n::T1,sz::T2) where {T1<:Real,T2<:Real}=new((Float64(n),Float64(sz)))
end

fieldnames(::Type{SPQN})=("N","Sz")
periods(::Type{SPQN})=(nothing,nothing)

"""
    N(qn:SPQN)

Get the particle number.
"""
N(qn::SPQN)=qn.values[1]

"""
    Sz(qn:SPQN)

Get the spin z-component.
"""
Sz(qn::SPQN)=qn.values[2]
