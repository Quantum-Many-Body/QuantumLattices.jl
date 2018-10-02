module GoodQuantumNumber
    import Base: ==,hash,show,length,zero,eltype,fieldnames,getproperty
    import Printf: @printf

    include("QuantumNumber.jl")
    include("QuantumNumbers.jl")
    include("QuantumNumberPack.jl")
end #module
