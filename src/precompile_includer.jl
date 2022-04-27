#! format: off
should_precompile =     true


# Don't edit the following! Instead change the script for `snoop_bot`.
ismultios = true
ismultiversion = true
# precompile_enclosure
@static if !should_precompile
    # nothing
elseif !ismultios && !ismultiversion
    @static if isfile(joinpath(@__DIR__, "../deps/SnoopCompile/precompile/precompile_QuantumLattices.jl"))
        include("../deps/SnoopCompile/precompile/precompile_QuantumLattices.jl")
        _precompile_()
    end
else
    @static if Sys.islinux() 
    @static if v"1.6.0-DEV" <= VERSION <= v"1.6.9" 
    @static if isfile(joinpath(@__DIR__, "../deps/SnoopCompile/precompile/linux/1.6/precompile_QuantumLattices.jl"))
        include("../deps/SnoopCompile/precompile/linux/1.6/precompile_QuantumLattices.jl")
        _precompile_()
    end
elseif v"1.7.0-DEV" <= VERSION <= v"1.7.9" 
    @static if isfile(joinpath(@__DIR__, "../deps/SnoopCompile/precompile/linux/1.7/precompile_QuantumLattices.jl"))
        include("../deps/SnoopCompile/precompile/linux/1.7/precompile_QuantumLattices.jl")
        _precompile_()
    end
else 
    end

else 
    end

end # precompile_enclosure
