module GoodQuantumNumber
    import Base: ==,hash,show
    import Printf: @printf
    export Qndtypes,QuantumNumber,SPQN,QuantumNumbers

    "The allowed data types for quantum numbers."
    const Qndtypes=Union{Int64,Float64}

    include("QuantumNumber.jl")
    include("QuantumNumbers.jl")
end #module
