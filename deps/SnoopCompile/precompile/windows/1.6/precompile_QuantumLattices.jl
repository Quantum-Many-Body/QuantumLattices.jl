# Use
#    @warnpcfail precompile(args...)
# if you want to be warned when a precompile directive fails
macro warnpcfail(ex::Expr)
    modl = __module__
    file = __source__.file === nothing ? "?" : String(__source__.file)
    line = __source__.line
    quote
        $(esc(ex)) || @warn """precompile directive
     $($(Expr(:quote, ex)))
 failed. Please report an issue in $($modl) (after checking for duplicates) or remove this directive.""" _file=$file _line=$line
    end
end


const __bodyfunction__ = Dict{Method,Any}()

# Find keyword "body functions" (the function that contains the body
# as written by the developer, called after all missing keyword-arguments
# have been assigned values), in a manner that doesn't depend on
# gensymmed names.
# `mnokw` is the method that gets called when you invoke it without
# supplying any keywords.
function __lookup_kwbody__(mnokw::Method)
    function getsym(arg)
        isa(arg, Symbol) && return arg
        @assert isa(arg, GlobalRef)
        return arg.name
    end

    f = get(__bodyfunction__, mnokw, nothing)
    if f === nothing
        fmod = mnokw.module
        # The lowered code for `mnokw` should look like
        #   %1 = mkw(kwvalues..., #self#, args...)
        #        return %1
        # where `mkw` is the name of the "active" keyword body-function.
        ast = Base.uncompressed_ast(mnokw)
        if isa(ast, Core.CodeInfo) && length(ast.code) >= 2
            callexpr = ast.code[end-1]
            if isa(callexpr, Expr) && callexpr.head == :call
                fsym = callexpr.args[1]
                if isa(fsym, Symbol)
                    f = getfield(fmod, fsym)
                elseif isa(fsym, GlobalRef)
                    if fsym.mod === Core && fsym.name === :_apply
                        f = getfield(mnokw.module, getsym(callexpr.args[2]))
                    elseif fsym.mod === Core && fsym.name === :_apply_iterate
                        f = getfield(mnokw.module, getsym(callexpr.args[3]))
                    else
                        f = getfield(fsym.mod, fsym.name)
                    end
                else
                    f = missing
                end
            else
                f = missing
            end
        else
            f = missing
        end
        __bodyfunction__[mnokw] = f
    end
    return f
end

function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{QuantumLattices.Essentials.DegreesOfFreedom.var"##s195#47",Any,Any,Any})
    Base.precompile(Tuple{QuantumLattices.Essentials.Frameworks.var"##s195#31",Any,Any,Any,Any,Any})
    Base.precompile(Tuple{QuantumLattices.Essentials.QuantumOperators.var"##s17#4",Any,Any,Any,Any,Any,Any})
    Base.precompile(Tuple{QuantumLattices.Essentials.QuantumOperators.var"##s17#7",Any,Any})
    Base.precompile(Tuple{QuantumLattices.Essentials.Spatials.var"##s164#33",Any,Any,Any,Type,Any})
    Base.precompile(Tuple{QuantumLattices.Essentials.Spatials.var"##s195#80",Any,Any,Any,Any})
    Base.precompile(Tuple{Type{Bonds},Lattice{1, PID, Float64}})
    Base.precompile(Tuple{Type{Generator},Tuple{Term{:Hopping, :t, Float64, Int64, TermCouplings{OperatorSum{Coupling{Int64, Tuple{FID{:*, Symbol, Symbol, Symbol}, FID{:*, Symbol, Symbol, Symbol}}, Subscripts{Tuple{NamedTuple{(:orbital, :spin), Tuple{Subscript{Tuple{Symbol, Symbol}}, Subscript{Tuple{Symbol, Symbol}}}}}, Tuple{Pair{UnitRange{Int64}, Tuple{String, String}}}}}, Tuple{CompositeIID{Tuple{FID{:*, Symbol, Symbol, Symbol}, FID{:*, Symbol, Symbol, Symbol}}}, Subscripts{Tuple{NamedTuple{(:orbital, :spin), Tuple{Subscript{Tuple{Symbol, Symbol}}, Subscript{Tuple{Symbol, Symbol}}}}}, Tuple{Pair{UnitRange{Int64}, Tuple{String, String}}}}}}, Nothing}, TermAmplitude{Nothing}, TermModulate{Val{false}, :t}}, Term{:Hubbard, :U, Float64, Int64, TermCouplings{OperatorSum{Coupling{Int64, NTuple{4, FID{:*, Symbol, Int64, Int64}}, Subscripts{Tuple{NamedTuple{(:orbital, :spin), Tuple{Subscript{NTuple{4, Symbol}}, Subscript{NTuple{4, Int64}}}}}, Tuple{Pair{UnitRange{Int64}, Tuple{String, String}}}}}, Tuple{CompositeIID{NTuple{4, FID{:*, Symbol, Int64, Int64}}}, Subscripts{Tuple{NamedTuple{(:orbital, :spin), Tuple{Subscript{NTuple{4, Symbol}}, Subscript{NTuple{4, Int64}}}}}, Tuple{Pair{UnitRange{Int64}, Tuple{String, String}}}}}}, Nothing}, TermAmplitude{Nothing}, TermModulate{Val{false}, :U}}},Bonds{(QuantumLattices.Essentials.Spatials.ZerothBonds(), QuantumLattices.Essentials.Spatials.InsideBonds(), QuantumLattices.Essentials.Spatials.AcrossBonds()), Lattice{1, PID, Float64}, Tuple{Vector{Point{1, PID, Float64}}, Vector{Bond{1, PID, Float64}}, Vector{Bond{1, PID, Float64}}}, AbstractBond{1, PID, Float64, R} where R},Hilbert{Fock{:f}, PID}})
    Base.precompile(Tuple{Type{Lattice},Symbol,Vector{Point{1, PID, Float64}}})
    Base.precompile(Tuple{typeof(QuantumLattices.Prerequisites.Traits._order),Val{:cid},Val{(:value, :cid, :subscripts)}})
    Base.precompile(Tuple{typeof(QuantumLattices.Prerequisites.Traits._order),Val{:id},Val{(:value, :id)}})
    Base.precompile(Tuple{typeof(QuantumLattices.Prerequisites.Traits._parametertype),Type{OperatorPack{Int64, Tuple{CompositeIID{Tuple{FID{:*, Symbol, Symbol, Symbol}, FID{:*, Symbol, Symbol, Symbol}}}, Subscripts{Tuple{NamedTuple{(:orbital, :spin), Tuple{Subscript{Tuple{Symbol, Symbol}}, Subscript{Tuple{Symbol, Symbol}}}}}, Tuple{Pair{UnitRange{Int64}, Tuple{String, String}}}}}}},Val{2}})
    Base.precompile(Tuple{typeof(QuantumLattices.Prerequisites.Traits._supertype),Type{Coupling{Int64, Tuple{FID{:*, Symbol, Symbol, Symbol}, FID{:*, Symbol, Symbol, Symbol}}, Subscripts{Tuple{NamedTuple{(:orbital, :spin), Tuple{Subscript{Tuple{Symbol, Symbol}}, Subscript{Tuple{Symbol, Symbol}}}}}, Tuple{Pair{UnitRange{Int64}, Tuple{String, String}}}}}},Val{:OperatorPack}})
    Base.precompile(Tuple{typeof(expand),Type{Lattice{1, PID, Float64}},Val{(QuantumLattices.Essentials.Spatials.AcrossBonds(),)}})
    Base.precompile(Tuple{typeof(expand),Type{Lattice{1, PID, Float64}},Val{(QuantumLattices.Essentials.Spatials.AllBonds(),)}})
    Base.precompile(Tuple{typeof(iterate),QuantumLattices.Prerequisites.SimpleTrees.TreeKeys{QuantumLattices.Prerequisites.SimpleTrees.SimpleTreeDepth, LatticeBonds, Nothing, SimpleTree{LatticeBonds, Nothing}},Vector{LatticeBonds}})
    Base.precompile(Tuple{typeof(iterate),QuantumLattices.Prerequisites.SimpleTrees.TreeKeys{QuantumLattices.Prerequisites.SimpleTrees.SimpleTreeDepth, LatticeBonds, Nothing, SimpleTree{LatticeBonds, Nothing}}})
    Base.precompile(Tuple{typeof(latticebondsstructure),Type{Lattice{1, PID, Float64}}})
    Base.precompile(Tuple{typeof(tile),Matrix{Float64},Vector{Vector{Float64}},Vector{Tuple{}}})
    let fbody = try __lookup_kwbody__(which(replace, (Term,))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Any,typeof(replace),Term,))
    end
end
end
