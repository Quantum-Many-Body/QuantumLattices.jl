var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": "CurrentModule=Hamiltonian"
},

{
    "location": "index.html#Hamiltonian-1",
    "page": "Home",
    "title": "Hamiltonian",
    "category": "section",
    "text": "Julia package for constructing and solving the Hamiltonians of quantum lattice systems.We provide a general framework to construct the symbolic representation of the Hamiltonian of any quantum lattice system, with the inputs as simple as its description by natural language. Based on this symbolic representation, we implement several algorithms, such as TBA, ED, CPT/VCA, DMRG, etc., to solve the quantum lattice system. Generic interfaces are offered to access to these algorithms. Only minor modifications need be made when the user alters from an algorithm to another."
},

{
    "location": "index.html#Package-Features-1",
    "page": "Home",
    "title": "Package Features",
    "category": "section",
    "text": "Unitcell Description Framework: by telling the information of the quantum lattice system within a unitcell, the construction of the symbolic representation of the Hamiltonian is just as simple as describing the system in a usual research paper;\nGeneric Engine-App Interfaces: by regarding the relation between algorithms and tasks as that between engines and apps, automatic project management is realized, including that of result recording, data caching, parameter updating, code logging, task dependency, etc, furthermore, all algorithms are initialized in quite similiar ways with only minor modifications needed."
},

{
    "location": "index.html#Supported-Algorithms-1",
    "page": "Home",
    "title": "Supported Algorithms",
    "category": "section",
    "text": "TBA: tight-binding approximation for fermionic/bosonic systems;\nED: exact diagonalizaiton for fermionic/hard-core-bosonic/spin systems;\nCPT/VCA: cluster perturbation theory and variational cluster approach for fermionic systems;\nDMRG: density matrix renormalization group for fermionic/hard-core-bosonic/spin systems;\nFBFM: spin wave theory for flatband ferromagnets."
},

{
    "location": "index.html#Python-counterpart-1",
    "page": "Home",
    "title": "Python counterpart",
    "category": "section",
    "text": "HamiltonianPy: in fact, the authors of this Julia package worked on the python package at first and only turned to Julia later."
},

{
    "location": "tutorial/UnitcellDescription.html#",
    "page": "Unitcell Description",
    "title": "Unitcell Description",
    "category": "page",
    "text": ""
},

{
    "location": "tutorial/UnitcellDescription.html#Unitcell-Description-1",
    "page": "Unitcell Description",
    "title": "Unitcell Description",
    "category": "section",
    "text": "This is the basic idea behind our framework, i.e. we only describe how the quantum lattice system looks like in a unitcell. To achieve this goal, we divide the problem into three steps:STEP 1: define the unitcell of the lattice;\nSTEP 2: define the local internal Hilbert spaces of the system on the unitcell;\nSTEP 3: define the terms connecting local Hilbert spaces on the sites of the unitcell."
},

{
    "location": "tutorial/EngineAppInterface.html#",
    "page": "Engine App Interface",
    "title": "Engine App Interface",
    "category": "page",
    "text": ""
},

{
    "location": "tutorial/EngineAppInterface.html#Engine-App-Interface-1",
    "page": "Engine App Interface",
    "title": "Engine App Interface",
    "category": "section",
    "text": "Althogh we can get the symbolic representation of the Hamiltonian by our unitcell-description framework, there still remains a long way to implement concrete algorithms such as TBA, ED, etc. Despite the quite different techinal details, algorithms shares common functionalities to be furnished with:provide tasks to be conducted with controlling parameters;\nrecord the results of some tasks for later use or analysis;\nupdate some parameters of the Hamiltonian to reconduct tasks;\nkeep logs during code executions for debug;\ncahe intermediate data to improve efficienty;\n...We thus provide a set of generic interfaces to resolve these problems, basiscally in the so called Engine-App mode. Specifically, algorithms are treated as Engine and taks as App. Engine deals with the cores of algorithms along with file mangament, parameter updating and data caching, while App decides the concrete tasks to be conducted and provides hyper controlling parameters."
},

{
    "location": "man/Utilities/Introduction.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Utilities"
},

{
    "location": "man/Utilities/Introduction.html#Introduction-1",
    "page": "Introduction",
    "title": "Introduction",
    "category": "section",
    "text": "This module contains the Utilities of the Hamiltonian package, all of whose variables will NOT be exported by Hamiltonian."
},

{
    "location": "man/Utilities/Introduction.html#Hamiltonian.Utilities.atol",
    "page": "Introduction",
    "title": "Hamiltonian.Utilities.atol",
    "category": "constant",
    "text": "Absolute tolerance for float numbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Introduction.html#Hamiltonian.Utilities.rtol",
    "page": "Introduction",
    "title": "Hamiltonian.Utilities.rtol",
    "category": "constant",
    "text": "Relative tolerance for float numbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Introduction.html#Hamiltonian.Utilities.Float",
    "page": "Introduction",
    "title": "Hamiltonian.Utilities.Float",
    "category": "type",
    "text": "Default float type.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Introduction.html#Hamiltonian.Utilities.FOrder",
    "page": "Introduction",
    "title": "Hamiltonian.Utilities.FOrder",
    "category": "type",
    "text": "Fortran order.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Introduction.html#Hamiltonian.Utilities.COrder",
    "page": "Introduction",
    "title": "Hamiltonian.Utilities.COrder",
    "category": "type",
    "text": "C order.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Introduction.html#Hamiltonian.Utilities.ind2sub",
    "page": "Introduction",
    "title": "Hamiltonian.Utilities.ind2sub",
    "category": "function",
    "text": "ind2sub(dims::Tuple,ind::Int,::Type{FOrder}) -> Tuple\nind2sub(dims::Tuple,ind::Int,::Type{COrder}) -> Tuple\n\nConvert an linear index to Cartesian index. Fortran-order or C-order can be assigned.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Introduction.html#Hamiltonian.Utilities.sub2ind",
    "page": "Introduction",
    "title": "Hamiltonian.Utilities.sub2ind",
    "category": "function",
    "text": "sub2ind(dims::NTuple{N,Int},inds::NTuple{N,Int},::Type{FOrder}) where N -> Int\nsub2ind(dims::NTuple{N,Int},inds::NTuple{N,Int},::Type{COrder}) where N -> Int\n\nConvert an Cartesian index to linear index. Fortran-order or C-order can be assigned.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Introduction.html#Hamiltonian.Utilities.decimaltostr",
    "page": "Introduction",
    "title": "Hamiltonian.Utilities.decimaltostr",
    "category": "function",
    "text": "decimaltostr(number::Integer,n::Int=5)\ndecimaltostr(number::Rational,n::Int=5)\ndecimaltostr(number::AbstractFloat,n::Int=5)\ndecimaltostr(number::Complex,n::Int=5)\n\nConvert a number to a string with at most n decimal places.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Introduction.html#Hamiltonian.Utilities.ordinal",
    "page": "Introduction",
    "title": "Hamiltonian.Utilities.ordinal",
    "category": "function",
    "text": "ordinal(number::Interger)\n\nConvert a positive number to its corresponding ordinal.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Introduction.html#Useful-constants-and-functions-1",
    "page": "Introduction",
    "title": "Useful constants and functions",
    "category": "section",
    "text": "atol\nrtol\nFloat\nFOrder\nCOrder\nind2sub\nsub2ind\ndecimaltostr\nordinal"
},

{
    "location": "man/Utilities/Introduction.html#Hamiltonian.Utilities.dimension",
    "page": "Introduction",
    "title": "Hamiltonian.Utilities.dimension",
    "category": "function",
    "text": "Generic interface of the dimension of some types.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Introduction.html#Generic-functions-for-overloading-1",
    "page": "Introduction",
    "title": "Generic functions for overloading",
    "category": "section",
    "text": "dimension"
},

{
    "location": "man/Utilities/Introduction.html#Prerequisites-for-Essentials-1",
    "page": "Introduction",
    "title": "Prerequisites for Essentials",
    "category": "section",
    "text": "Pages=[\n    \"Factory.md\",\n    \"Tree.md\",\n    \"NamedVector.md\",\n    ]\nDepth=2"
},

{
    "location": "man/Utilities/Introduction.html#Necessities-for-Algorithms-1",
    "page": "Introduction",
    "title": "Necessities for Algorithms",
    "category": "section",
    "text": "Pages=[\n    \"QuantumNumber.md\",\n    ]\nDepth=2"
},

{
    "location": "man/Utilities/Factory.html#",
    "page": "Factory",
    "title": "Factory",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Utilities.Factorypush!(LOAD_PATH,\"../../../../src/\")\nusing Hamiltonian.Utilities.Factory"
},

{
    "location": "man/Utilities/Factory.html#Factory-1",
    "page": "Factory",
    "title": "Factory",
    "category": "section",
    "text": "The aim of Factory is to provide tools to hack into Julia codes without knowing the details of their abstract syntax trees and regularize the mechanism to \"escape\" variables in Expr expressions, so that users can manipulate the existing codes, modify them and generate new ones in macros. In particular, Factory in this module means the representation of certain blocks of Julia codes by a usual Julia struct. This representation is much easier to comprehend than the canonical Expr representation and makes it far more convenient to define macros. In general, we propose the following requirements that any factory must satisfy:DECOMPOSITION - An Expr expression can be decomposed into its corresponding factory by the factory\'s constructor.\nCOMPOSITION - A factory can compose its corresponding Expr expression by calling itself.\nESCAPE - A variable should be escaped or not in the composed Expr expression by a factory depends on keyword arguments escaped or unescaped passed to the factory call, which are tuple of Symbols. Specifically, when escaped is provided, a variable should be escaped if its name is in escaped, whereas, when unescaped is provided, a variable should be escaped if its name is not in unescaped. A factory can choose which keyword arguments should be used for convenience, or can choose both with different keyword arguments for different parts of it.The first two requirements define the basic interfaces to interact with factories, and the third requirement proposes the escape mechanism of variables.Out of practical purposes, we implemente 7 kinds of factories, i.e. Inference, Argument, Parameter, Field, Block, FunctionFactory and TypeFactory, which represent a type inference, a function argument, a method or type parameter, a struct field, a begin ... end block, a function itself and a struct itself, respectively. Some of the basic methods making the above requirements fulfilled with these types are based on the powerful functions defined in MacroTools."
},

{
    "location": "man/Utilities/Factory.html#Inference-1",
    "page": "Factory",
    "title": "Inference",
    "category": "section",
    "text": "An Inference has 3 attributes:head::Union{Symbol,Nothing}: the head of the type inference, which must be one of (nothing,:(::),:(<:),:curly)\nname::Union{Symbol,Nothing}: the name of the type inference\nparams::Union{Inference,Vector{Inference},Nothing}: the parameters of the type inferenceAll valid expressions representing type inferences can be passed to the constructor:Inference(:T)\nInference(:(::T))\nInference(:(<:Number))\nInference(:(Vector{T}))\nInference(:(Vector{Tuple{String,Int}}))\nInference(:(Type{<:Number}))On the other hand, you can use the macro @inference to construct an Inference directly from a type inference:@inference Vector{Tuple{String,Int}}note: Note\nInference is a recursive struct, i.e. it recursively decomposes a type inference until the final type inference is just a Symbol.\nWhen the input expression is a Symbol, the head and params attributes of the resulting Inference is nothing. Otherwise, its head is the same with that of the input expression, and the args of the input expression will be further decomposed, whose result will be stored in params.\nWhen the head of the input expression is :(::) or :(<:), the params is an Inference whereas when the head of the input expression is :curly, the params is a Vector{Inference}.Inference uses the keyword argument unescaped to escape variables, e.g.Inference(:(Vector{T}))() |> println\nInference(:(Vector{T}))(unescaped=(:T,)) |> println\nInference(:(Vector{T}))(unescaped=(:Vector,:T)) |> println"
},

{
    "location": "man/Utilities/Factory.html#Argument-1",
    "page": "Factory",
    "title": "Argument",
    "category": "section",
    "text": "An Argument has 4 attributes:name::Union{Symbol,Nothing}: the name of the argument\ntype::Inference: the type inference of the argument\nslurp::Bool: whether the argument should be expanded by ...\ndefault::Any: the default value of the argument, nothing for those with no default valuesAll valid expressions representing the arguments of functions can be passed to the constructor:Argument(:arg)\nArgument(:(arg::ArgType))\nArgument(:(arg::ArgType...))\nArgument(:(arg::ArgType=default))Or you can use the macro @argument for a direct construction from an argument declaration:@argument arg::ArgType=defaultThe construction from such expressions is based on the the MacroTools.splitarg function.Argument also uses the keyword argument unescaped to escape variables, e.g.Argument(:(arg::ArgType=default))(unescaped=(:ArgType,:default)) |> printlnIt can be seen the name of an argument will never be escaped even though it is not in the unescaped tuple. This is obvious since the name of a function argument is always local. By the way, the composition of an Argument expression is based on the MacroTools.combinearg function."
},

{
    "location": "man/Utilities/Factory.html#Parameter-1",
    "page": "Factory",
    "title": "Parameter",
    "category": "section",
    "text": "A Parameter has 2 attributes:name::Union{Symbol,Nothing}: the name of the parameter\ntype::Union{Inference,Nothing}: the type inference of the parameterAll expressions that represent type parameters or method parameters are allowed to be passed to the constructor:Parameter(:T)\nParameter(:(<:Number))\nParameter(:(T<:Number))\nParameter(:(::Int))The macro @parameter completes the construction directly from a parameter declaration:@parameter T<:Numbernote: Note\nWe use nothing to denote a missing name or type.\nTwo subtle situations of type/method parameters, e.g. MyType{T} and MyType{Int}, are distinguished by Parameter(:T) and Parameter(:(::Int)). The name and type attributes of the resulting Parameters are, for the first case, :T and nothing, while for the second case, nothing and :T, respectively. Moreover, the callings of the factories for these two cases are also different, e.g. Parameter(:T)()==:T and Parameter(:(::Int))==:(<:$(esc(Int))).Parameter uses the keyword argument unescaped to escape variables, too, e.g.Parameter(:(N<:Vector{T}))(unescaped=(:T,)) |> printlnAs is similar to Argument, the name of a method/type parameter will never be escaped because of its local scope."
},

{
    "location": "man/Utilities/Factory.html#Field-1",
    "page": "Factory",
    "title": "Field",
    "category": "section",
    "text": "A Field has 2 attributes:name::Symbol: the name of the field\ntype::Inference: the type inference of the fieldLegal expressions can be used to construct a Field instance by its constructor:Field(:field)\nField(:(field::FieldType))\nField(:(field::ParametricType{T}))The macro @field is also provided to help the construction directly from a field declaration:@field field::FieldTypeThe construction from these expressions is based on the MacroTools.splitarg function.Field uses the keyword argument unescaped to escape variables as well, e.g.Field(:(field::Dict{N,D}))(unescaped=(:N,:D)) |> printlnThe name of a struct will never be escaped either because it is a local variable tightly binding to a struct. It is noted that the composition of field expressions is based on the MacroTools.combinefield function."
},

{
    "location": "man/Utilities/Factory.html#Block-1",
    "page": "Factory",
    "title": "Block",
    "category": "section",
    "text": "A Block has only one attribute:body::Vector{Any}: the body of the begin ... end blockAny expression can be passed to the constructor of Block:Block(:(x=1))\nBlock(:(x=1;y=2))\nBlock(:(begin x=1 end))\nBlock(quote\n        x=1\n        y=2\n    end)Or you can construct a Block instance directly from any code by the macro @block:@block x=1 y=2The body of a block can also be extended by the push! function or the @push! macro.note: Note\nThe body of a Block is somewhat \"flattened\", i.e. it contains no begin ... end blocks. During the initialization, any such input block will be unblocked and added to the body part by part. So is the push! and @push! procedures.\nAll LineNumberNodes generated by the input codes will also be included in the block\'s body. However, you can use rmlines! or @rmlines! to remove them from the body of an existing Block, or use rmlines or @rmlines to get a copy with them removed in the body.Different from previous factories, Block uses the keyword argument escaped to escape variables. This is because variables in a block are often local ones and should not be escaped. Therefore, only those defined in other modules should be noted and escaped, which usually constitute the minority. For example,Block(:(x=1;y=2;z=Int[1,2,3]))(escaped=(:Int,)) |> println"
},

{
    "location": "man/Utilities/Factory.html#FunctionFactory-1",
    "page": "Factory",
    "title": "FunctionFactory",
    "category": "section",
    "text": "A FunctionFactory has 7 attributes:name::Union{Symbol,Expr}: the name of the function\nparams::Vector{Inference}: the method parameters of the function\nargs::Vector{Argument}: the positional arguments of the function\nkwargs::Vector{Argument}: the keyword arguments of the function\nrtype::Inference: the return type of the function\nwhereparams::Vector{Parameter}: the method parameters specified by the where keyword\nbody::Block: the body of the functionAll expressions that represent functions are allowed to be passed to the constructor:FunctionFactory(:(f()=nothing))\nFunctionFactory(:(f(x)=x))\nFunctionFactory(:(f(x::Int,y::Int;choice::Function=sum)=choice(x,y)))\nFunctionFactory(:(f(x::T,y::T;choice::Function=sum) where T<:Number=choice(x,y)))\nFunctionFactory(:((f(x::T,y::T;choice::Function=sum)::T) where T<:Number=choice(x,y)))\nFunctionFactory(:(\n    function (f(x::T,y::T;choice::Function=sum)::T) where T<:Number\n        choice(x,y)\n    end\n))\nFunctionFactory(\n    quote\n        function (f(x::T,y::T;choice::Function=sum)::T) where T<:Number\n            choice(x,y)\n        end\n    end\n)Similarly, an instance can also be constructed from the macro @functionfactory:@functionfactory (f(x::T,y::T;choice::Function=sum)::T) where T<:Number=choice(x,y)The construction from such expressions are based on the MacroTools.splitdef function.note: Note\nSince Julia 0.7, the form MyType{D}(data::D) where D only appears in struct constructors, therefore, the attribute :params of a function factory is nonempty only when this factory aims to represent a struct constructor.\nUsually, the name of a function factory is a Symbol. However, if the factory aims to extend some methods of a function defined in another module, e.g., Base.eltype, the name will be an Expr.Since FunctionFactory involves not only factories using unescaped but also factories using escaped, it adopts both to escape variables, with unescaped for params, args, kwargs, rtype and whereparams while escaped for name and body. It is worth to emphasize that the name of a function factory is affected by the escaped argument. Specifically, when the name is a Symbol and is in the escaped tuple, it will be escaped. Otherwise it will not, especially when it is an Expr, it will never be escaped because an Expr cannot be a element of a NTuple{N,Symbol} where N. See examples,FunctionFactory(:((f(x::T,y::T;choice::Function=sum)::T) where T<:Number=max(x,y,choice(x,y))))(unescaped=(:T,),escaped=(:f,:max,)) |> println\nFunctionFactory(:((f(x::T,y::T;choice::Function=sum)::T) where T<:Number=max(x,y,choice(x,y))))(unescaped=(:T,),escaped=(:max,)) |> printlnThe compositions of function expressions are based on the MacroTools.combinedef function.Other features include:Positional arguments can be added by addargs! or @addargs!\nKeyword arguments can be added by addkwargs! or @addkwargs!\nWhere parameters can be added by addwhereparams! or @addwhereparams!\nBody can be extended by extendbody! or @extendbody!"
},

{
    "location": "man/Utilities/Factory.html#TypeFactory-1",
    "page": "Factory",
    "title": "TypeFactory",
    "category": "section",
    "text": "A TypeFactory has 6 attributes:name::Symbol: the name of the struct\nmutable::Bool: whether or not the struct is mutable\nparams::Vector{Parameter}: the type parameters of the struct\nsupertype::Inference: the supertype of the struct\nfields::Vector{Field}: the fields of the struct\nconstructors::Vector{FunctionFactory}: the inner constructors of the structAny expression representing valid struct definitions can be passed to the constructor:TypeFactory(:(struct StructName end))\nTypeFactory(:(struct StructName{T} end))\nTypeFactory(:(struct Child{T} <: Parent{T} end))\nTypeFactory(:(\n    struct Child{T<:Number} <: Parent{T}\n        field1::T\n        field2::T\n    end\n))\nTypeFactory(\n    quote\n        struct Child{T<:Number} <: Parent{T}\n            field1::T\n            field2::T\n        end\n    end\n)Also, the macro @typefactory supports the construction directly from a type definition:@typefactory struct Child{T<:Number} <: Parent{T}\n                field1::T\n                field2::T\n                Child(field1::T,field2::T=zero(T)) where T=new{T}(field1,field2)\n            endThe construction from these expressions is based on the MacroTools.splitstructdef function.TypeFactory also uses both keyword arguments, i.e. unescaped and escaped, to escape variables, with the former for params, supertype and fields, the latter for name, and both for constructors. For example,(@typefactory struct Child{T<:Number} <: Parent{T} field::T; Child(field::T) where T=new{T}(field) end)(unescaped=(:T,),escaped=(:Child,)) |>printlnThe composition of a type expression is based on the MacroTools.combinestructdef function.Other features include:Fields can be added by addfields! or @addfields!\nType parameters can be added by addparams! or @addparams!\nInner constructors can be added by addconstructors! or @addconstructors!"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.FExpr",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.FExpr",
    "category": "constant",
    "text": "Factory expression types, which is defined as Union{Symbol,Expr}.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.RawExpr",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.RawExpr",
    "category": "constant",
    "text": "Whether or not to show raw expressions of factoies.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.Argument",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.Argument",
    "category": "type",
    "text": "Argument(name::Union{Symbol,Nothing},type::Inference,slurp::Bool,default::Any)\nArgument(;name::Union{Symbol,Nothing}=nothing,type::Inference=Inference(:Any),slurp::Bool=false,default::Any=nothing)\nArgument(expr::FExpr)\n\nThe struct to describe a argument of a function.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.Argument-Tuple{Val{true}}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.Argument",
    "category": "method",
    "text": "(a::Argument)(::typeof(RawExpr)) -> Expr\n(a::Argument)(;unescaped::NTuple{N,Symbol}=()) where N -> Expr\n\nConvert an Argument to the Expr representation of the argument it describes.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.Block",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.Block",
    "category": "type",
    "text": "Block(parts::FExpr...)\n\nThe struct to describe a begin ... end block.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.Block-Tuple{Val{true}}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.Block",
    "category": "method",
    "text": "(b::Block)(::typeof(RawExpr)) -> Expr\n(b::Block)(;escaped::NTuple{N,Symbol}=()) where N -> Expr\n\nConvert a Block to the Expr representation of the begin ... end block it describes.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.Field",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.Field",
    "category": "type",
    "text": "Field(name::Symbol,type::Inference)\nField(;name::Symbol,type::FExpr=Inference(:Any))\nField(expr::FExpr)\n\nThe struct to describe a field of a struct.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.Field-Tuple{Val{true}}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.Field",
    "category": "method",
    "text": "(f::Field)(::typeof(RawExpr)) -> Expr\n(f::Field)(;unescaped::NTuple{N,Symbol}=()) where N -> Expr\n\nConvert a Field to the Expr representation of the field it describes.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.FunctionFactory",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.FunctionFactory",
    "category": "type",
    "text": "FunctionFactory(name::FExpr,params::Vector{Inference},args::Vector{Argument},kwargs::Vector{Argument},rtype::Inference,whereparams::Vector{Parameter},body::Block)\nFunctionFactory(    ;name::FExpr,\n                    params::Vector{Inference}=Inference[],\n                    args::Vector{Argument}=Argument[],\n                    kwargs::Vector{Argument}=Argument[],\n                    rtype::Inference=Inference(:Any),\n                    whereparams::Vector{Parameter}=Parameter[],\n                    body::Block=Block()\n                    )\nFunctionFactory(expr::Expr)\n\nThe struct to describe a function.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.FunctionFactory-Tuple{Val{true}}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.FunctionFactory",
    "category": "method",
    "text": "(ff::FunctionFactory)(::typeof(RawExpr)) -> Expr\n(ff::FunctionFactory)(;unescaped::NTuple{N,Symbol}=(),escaped::NTuple{M,Symbol}=()) where {N,M} -> Expr\n\nConvert a FunctionFactory to the Expr representation of the function it describes.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.Inference",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.Inference",
    "category": "type",
    "text": "Inference(head::Union{Symbol,Nothing},name::Union{Symbol,Nothing},params::Union{Inference,Vector{Inference},Nothing})\nInference(;\n        head::Union{Symbol,Nothing}=nothing,\n        name::Union{Symbol,Nothing}=nothing,\n        params::Union{Inference,Vector{Inference},Nothing}=nothing,\n        )\nInference(expr::FExpr)\n\nThe struct to describe a type inference.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.Inference-Tuple{Val{true}}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.Inference",
    "category": "method",
    "text": "(i::Inference)(::typeof(RawExpr)) -> FExpr\n(i::Inference)(;unescaped::NTuple{N,Symbol}=()) where N -> FExpr\n\nConvert a Inference to the Expr representation of the type inference it describes.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.Parameter",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.Parameter",
    "category": "type",
    "text": "Parameter(name::Union{Symbol,Nothing},type::Union{Inference,Nothing})\nParameter(;name::Union{Symbol,Nothing}=nothing,type::Union{Inference,Nothing}=nothing)\nParameter(expr::FExpr)\n\nThe struct to describe a parameter of a function or a type.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.Parameter-Tuple{Val{true}}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.Parameter",
    "category": "method",
    "text": "(p::Parameter)(::typeof(RawExpr)) -> FExpr\n(p::Parameter)(;unescaped::NTuple{N,Symbol}=()) where N -> FExpr\n\nConvert a Parameter to the Expr representation of the parameter it describes.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.TypeFactory",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.TypeFactory",
    "category": "type",
    "text": "TypeFactory(name::Symbol,mutable::Bool,params::Vector{Parameter},supertype::Inference,fields::Vector{Field},constructors::Vector{FunctionFactory})\nTypeFactory(    ;name::Symbol,\n                mutable::Bool=false,\n                params::Vector{Parameter}=Parameter[],\n                supertype::Inference=Inference(:Any),\n                fields::Vector{Field}=Field[],\n                constructors::Vector{FunctionFactory}=FunctionFactory[],\n                )\nTypeFactory(expr::Expr)\n\nThe struct to describe a struct.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.TypeFactory-Tuple{Val{true}}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.TypeFactory",
    "category": "method",
    "text": "(tf::TypeFactory)(::typeof(RawExpr)) -> Expr\n(tf::TypeFactory)(;unescaped::NTuple{N,Symbol}=()) where N -> Expr\n\nConvert a TypeFactory to the Expr representation of the struct it describes.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.@addargs!-Tuple{Any,Vararg{Union{Expr, Symbol},N} where N}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.@addargs!",
    "category": "macro",
    "text": "@addargs! ff args::FExpr...\n\nAdd a couple of positional arguments to a function factory.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.@addconstructors!-Tuple{Any,Vararg{Expr,N} where N}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.@addconstructors!",
    "category": "macro",
    "text": "@addconstructors! tf constructors::Expr...\n\nAdd a couple of constructors to a type factory.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.@addfields!-Tuple{Any,Vararg{Union{Expr, Symbol},N} where N}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.@addfields!",
    "category": "macro",
    "text": "@addfields! tf fields::FExpr...\n\nAdd a couple of fields to a type factory.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.@addkwargs!-Tuple{Any,Vararg{Union{Expr, Symbol},N} where N}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.@addkwargs!",
    "category": "macro",
    "text": "@addkwargs! ff kwargs::FExpr...\n\nAdd a couple of keyword arguments to a function factory.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.@addparams!-Tuple{Any,Vararg{Union{Expr, Symbol},N} where N}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.@addparams!",
    "category": "macro",
    "text": "@addparams! f params::FExpr...\n\nAdd a couple of method parameters to a function factory or a type factory.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.@addwhereparams!-Tuple{Any,Vararg{Union{Expr, Symbol},N} where N}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.@addwhereparams!",
    "category": "macro",
    "text": "@addwhereparams! ff whereparams::FExpr...\n\nAdd a couple of method parameters to a function factory or a type factory.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.@argument-Tuple{Union{Expr, Symbol}}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.@argument",
    "category": "macro",
    "text": "@argument expr::FExpr\n\nConstruct an Argument directly from an argument statement.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.@block-Tuple{Vararg{Union{Expr, Symbol},N} where N}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.@block",
    "category": "macro",
    "text": "@block parts::FExpr...\n\nConstruct a Block directly from a begin ... end block definition.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.@extendbody!-Tuple{Any,Vararg{Union{Expr, Symbol},N} where N}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.@extendbody!",
    "category": "macro",
    "text": "@extendbody! ff parts::FExpr...\n\nExtend the body of a function factory.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.@field-Tuple{Union{Expr, Symbol}}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.@field",
    "category": "macro",
    "text": "@field expr::FExpr\n\nConstruct a Field directly from a field statement.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.@functionfactory-Tuple{Expr}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.@functionfactory",
    "category": "macro",
    "text": "@functionfactory expr::FExpr\n\nConstruct a FunctionFactory directly from a function definition.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.@inference-Tuple{Union{Expr, Symbol}}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.@inference",
    "category": "macro",
    "text": "@inference expr::FExpr\n\nConstruct an Inference directly from a type inference.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.@parameter-Tuple{Union{Expr, Symbol}}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.@parameter",
    "category": "macro",
    "text": "@parameter expr::FExpr\n\nConstruct a Parameter directly from an parameter statement.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.@push!-Tuple{Any,Vararg{Union{Expr, Symbol},N} where N}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.@push!",
    "category": "macro",
    "text": "@push! b parts::FExpr...\n\nPush other parts into the body of a block.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.@rmlines!-Tuple{Expr}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.@rmlines!",
    "category": "macro",
    "text": "@rmlines! b::Expr\n\nRemove line number nodes in the body of a block.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.@rmlines-Tuple{Expr}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.@rmlines",
    "category": "macro",
    "text": "@rmlines b::Expr\n\nReturn a copy of a block with the line number nodes removed.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.@typefactory-Tuple{Expr}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.@typefactory",
    "category": "macro",
    "text": "@typefactory expr::Expr\n\nConstruct a TypeFactory directly from a type definition.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.addargs!-Tuple{Hamiltonian.Utilities.Factory.FunctionFactory}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.addargs!",
    "category": "method",
    "text": "addargs!(ff::FunctionFactory,args::Argument...) -> FunctionFactory\naddargs!(ff::FunctionFactory,args::FExpr...) -> FunctionFactory\n\nAdd a couple of positional arguments to a function factory.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.addconstructors!-Tuple{Hamiltonian.Utilities.Factory.TypeFactory}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.addconstructors!",
    "category": "method",
    "text": "addconstructors!(tf::TypeFactory,constructors::FunctionFactory...) -> TypeFactory\naddconstructors!(tf::TypeFactory,constructors::Expr...) -> TypeFactory\n\nAdd a couple of constructors to a type factory.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.addfields!-Tuple{Hamiltonian.Utilities.Factory.TypeFactory}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.addfields!",
    "category": "method",
    "text": "addfields!(tf::TypeFactory,fields::Field...) -> TypeFactory\naddfields!(tf::TypeFactory,fields::FExpr...) -> TypeFactory\n\nAdd a couple of fields to a type factory.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.addkwargs!-Tuple{Hamiltonian.Utilities.Factory.FunctionFactory}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.addkwargs!",
    "category": "method",
    "text": "addkwargs!(ff::FunctionFactory,kwargs::Argument...) -> FunctionFactory\naddkwargs!(ff::FunctionFactory,kwargs::FExpr...) -> FunctionFactory\n\nAdd a couple of keyword arguments to a function factory.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.addparams!-Tuple{Union{FunctionFactory, TypeFactory}}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.addparams!",
    "category": "method",
    "text": "addparams!(f::FunctionFactory,params::Inference...) ->FunctionFactory\naddparams!(f::FunctionFactory,params::FExpr...) -> FunctionFactory\naddparams!(f::TypeFactory,params::Parameter...) -> TypeFactory\naddparams!(f::TypeFactory,params::FExpr...) -> TypeFactory\n\nAdd a couple of parameters to a function factory or a type factory.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.addwhereparams!-Tuple{Hamiltonian.Utilities.Factory.FunctionFactory}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.addwhereparams!",
    "category": "method",
    "text": "addwhereparams!(ff::FunctionFactory,whereparams::Parameter...) -> FunctionFactory\naddwhereparams!(ff::FunctionFactory,whereparams::FExpr...) -> FunctionFactory\n\nAdd a couple of method where parameters to a function factory or a type factory.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.escape-Union{Tuple{Any}, Tuple{N}, Tuple{Any,Tuple{Vararg{Symbol,N}}}} where N",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.escape",
    "category": "method",
    "text": "escape(expr,::NTuple{N,Symbol}=()) where N -> Any\nescape(expr::Symbol,escaped::NTuple{N,Symbol}=()) where N -> FExpr\nescape(expr::Expr,escaped::NTuple{N,Symbol}=()) where N -> Expr\n\nEscape the variables sepecified by escaped in the input expression.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.extendbody!-Tuple{Hamiltonian.Utilities.Factory.FunctionFactory}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.extendbody!",
    "category": "method",
    "text": "extendbody!(ff::FunctionFactory,parts::FExpr...) -> FunctionFactory\nextendbody!(ff::FunctionFactory,parts::Block...) -> FunctionFactory\n\nExtend the body of a function factory.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.rmlines!-Tuple{Hamiltonian.Utilities.Factory.Block}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.rmlines!",
    "category": "method",
    "text": "rmlines!(b::Block) -> Block\n\nRemove line number nodes in the body of a block.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#MacroTools.rmlines-Tuple{Hamiltonian.Utilities.Factory.Block}",
    "page": "Factory",
    "title": "MacroTools.rmlines",
    "category": "method",
    "text": "rmlines(b::Block) -> Block\n\nReturn a copy of a block with the line number nodes removed.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.AbstractFactory",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.AbstractFactory",
    "category": "type",
    "text": "AbstractFactory\n\nAbstract type for all concrete factories.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Base.:==-Union{Tuple{F}, Tuple{F,F}} where F<:Hamiltonian.Utilities.Factory.AbstractFactory",
    "page": "Factory",
    "title": "Base.:==",
    "category": "method",
    "text": "==(f1::F,f2::F) where F<:AbstractFactory -> Bool\n\nOverloaded == operator.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Base.names-Tuple{Any}",
    "page": "Factory",
    "title": "Base.names",
    "category": "method",
    "text": "names(expr::Union{Symbol,Expr}) -> Vector{Symbol}\n\nGet all the variable names in an Expr.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Base.push!-Tuple{Hamiltonian.Utilities.Factory.Block}",
    "page": "Factory",
    "title": "Base.push!",
    "category": "method",
    "text": "push!(b::Block,parts::FExpr...) -> Block\npush!(b::Block,parts::Block...) -> Block\n\nPush other parts into the body of a block.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Base.replace-Tuple{Hamiltonian.Utilities.Factory.AbstractFactory}",
    "page": "Factory",
    "title": "Base.replace",
    "category": "method",
    "text": "replace(f::AbstractFactory;kwargs...) -> typeof(f)\n\nReturn a copy of a concrete AbstractFactory with some of the field values replaced by the keyword arguments.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Base.show-Tuple{IO,Hamiltonian.Utilities.Factory.AbstractFactory}",
    "page": "Factory",
    "title": "Base.show",
    "category": "method",
    "text": "show(io::IO,f::AbstractFactory)\n\nShow a concrete AbstractFactory.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Manual-1",
    "page": "Factory",
    "title": "Manual",
    "category": "section",
    "text": "Modules=[Factory]\nOrder=  [:module,:constant,:type,:macro,:function]"
},

{
    "location": "man/Utilities/Tree.html#",
    "page": "Tree",
    "title": "Tree",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Utilities.Treepush!(LOAD_PATH,\"../../../../src/\")\nusing Hamiltonian.Utilities.Tree"
},

{
    "location": "man/Utilities/Tree.html#Tree-1",
    "page": "Tree",
    "title": "Tree",
    "category": "section",
    "text": "The aim of the tree in this module is to represent the standard tree structure in efficiency-non-sensitive cases. Please note that the default implementation of tree methods are far from optimal in efficiency. Therefore, please DO NOT use it if you need an efficient tree for addition, deletion, sort and inquiry. This module of codes apply only when the structure of tree matters but not the efficiency."
},

{
    "location": "man/Utilities/Tree.html#AbstractTree-1",
    "page": "Tree",
    "title": "AbstractTree",
    "category": "section",
    "text": "AbstractTree{N,D} is the abstract type for all concrete trees. By design, it has two type parameters:N: the type of the tree\'s node\nD: the type of the tree\'s dataTo fully utilize the methods designed for a tree structure, in our protocol, a concrete subtype must implement the following methods:inquiry related methods\neltype(tree::AbstractTree) -> NTuple\nroot(tree::AbstractTree{N,D}) where {N,D} -> Union{N,Nothing}\nhaskey(tree::AbstractTree{N,D},node::N) where {N,D} -> Bool\nlength(tree::AbstractTree) -> Int\nparent(tree::AbstractTree{N,D},node::N,superparent::Union{N,Nothing}=nothing) where {N,D} -> Union{N,Nothing}\nchildren(tree::AbstractTree{N,D},node::N) where {N,D} -> Vector{N}\nGet a tree\'s type parameters.\nGet a tree\'s root node (nothing for empty trees)\nGet the number of a tree\'s nodes.\nCheck whether a node is in a tree.\nGet the parent of a tree\'s node or return superparent when the input node is the tree\'s root.\nGet the children of a tree\'s node.\nstructure modification related methods\naddnode!(tree::AbstractTree{N,D},parent::Union{N,Nothing},node::N) where {N,D}\ndeletenode!(tree::AbstractTree{N,D},node::N) where {N,D}\nUpdate the structure of a tree by adding a node. When the parent is nothing, the input tree must be empty and the input node becomes the tree\'s root.\nUpdate the structure of a tree by deleting a node.\nindex related methods\ngetindex(tree::AbstractTree{N,D},node::N) where {N,D} -> D\nsetindex!(tree::AbstractTree{N,D},node::N,data::D) where {N,D}\nGet the data of a tree\'s node\nSet the data of a tree\'s node.Based on these methods, we implement several generic functions for inquiries and manipulationsinquiry for type parameters: keytype, valtype\nexpansion over nodes/data-records: keys, values, pairs\ninquiry for info of nodes: isleaf, level\ninquiry for nodes: ancestor, descendants, siblings, leaves\nmodification: push!, append!, delete!, empty!And optionally, when a subtype implement the following method,empty(tree::AbstractTree) -> typeof(tree)which constructs an empty tree of the same type with the input one, two more more methods are supported:subtree: Get a subtree starting from a node.\nmove!: Move a subtree to a new position."
},

{
    "location": "man/Utilities/Tree.html#TreeCore-and-SimpleTree-1",
    "page": "Tree",
    "title": "TreeCore and SimpleTree",
    "category": "section",
    "text": "To implement all the prerequisites listed above costs a bit efforts. We provide two lazy ways to get over this:Inheritance with a specific attribute TREECORE::TreeCore\nInclusion an attribute which is an instance of SimpleTree"
},

{
    "location": "man/Utilities/Tree.html#TreeCore-1",
    "page": "Tree",
    "title": "TreeCore",
    "category": "section",
    "text": "TreeCore{N,D}, as the literal meaning indicates, is the core of a tree. It encapsulates all the data structures needed by the default implementation, which constains 4 attributes:root::N: the tree\'s root node\ncontents::Dict{N,D}: the tree\'s (node,data) pairs\nparent::Dict{N,N}: records of the parent of each of the tree\'s nodes\nchildren::Dict{N,Vector{N}}: records of the children of each of the tree\'s nodesAs above, the first lazy way is to include this struct with the special attribute name :TREECORE in your concrete subtype. This process can be even lazier, in that we provide a macro @tree to decorate your \"raw\" struct automatically, e.g.@tree struct SimpleSubTree end\n@tree struct SubTreeWithTreeParameters end {N<:AbstractString,D<:Number}\n@tree struct SubTreeWithCertainTreeParameters end {::String,::Int}\n@tree struct SubTreeWithFields info::Vector{Int} end {N<:AbstractString,D<:Number}\n@tree struct SubTreeWithParametricFields{T} info::Vector{T} end {N<:AbstractString,D<:Number}\n@tree struct SubTreeWithOverlappedParametricFields{N} info::Vector{N} end {N<:AbstractString,D<:Number}"
},

{
    "location": "man/Utilities/Tree.html#SimpleTree-1",
    "page": "Tree",
    "title": "SimpleTree",
    "category": "section",
    "text": "SimpleTree{N,D} is the minimum struct that implements all the default tree methods. You can include an instance of it as an attribute in your own type to utilize all the tree methods."
},

{
    "location": "man/Utilities/Tree.html#Hamiltonian.Utilities.Tree.treedepth",
    "page": "Tree",
    "title": "Hamiltonian.Utilities.Tree.treedepth",
    "category": "constant",
    "text": "Depth first search or iteration.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Hamiltonian.Utilities.Tree.treewidth",
    "page": "Tree",
    "title": "Hamiltonian.Utilities.Tree.treewidth",
    "category": "constant",
    "text": "Width first search or iteration.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Hamiltonian.Utilities.Tree.AbstractTree",
    "page": "Tree",
    "title": "Hamiltonian.Utilities.Tree.AbstractTree",
    "category": "type",
    "text": "AbstractTree{Node,Data}\n\nAbstract type for all concrete trees.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Hamiltonian.Utilities.Tree.SimpleTree",
    "page": "Tree",
    "title": "Hamiltonian.Utilities.Tree.SimpleTree",
    "category": "type",
    "text": "SimpleTree()\n\nThe minimum tree structure that implements all the default tree methods.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Hamiltonian.Utilities.Tree.TreeCore",
    "page": "Tree",
    "title": "Hamiltonian.Utilities.Tree.TreeCore",
    "category": "type",
    "text": "TreeCore()\n\nThe core of a tree.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Hamiltonian.Utilities.Tree.@tree",
    "page": "Tree",
    "title": "Hamiltonian.Utilities.Tree.@tree",
    "category": "macro",
    "text": "@tree structdef treeparams::Union{Expr,Nothing}=nothing\n\nDecorate a \"raw\" struct to be a subtype of AbstractTree.\n\nnote: Note\nA \"raw\" struct means:\nIt has no explicit supertype;\nIt has no inner constructor;\nIt has no attribute :TREECORE.\nThe keytype and valtype can be assigned by the argument treeparams in the form {keytype,valtype}.\nWhen the formal argument names of keytype and valtype are not assigned, those of the type parameters of the raw struct, if any, cannot be :N or :D. For example, all of the following codes\n@tree struct SubTreeWithWrongTypeParameterNames{N} info::Vector{N} end\n@tree struct SubTreeWithWrongTypeParameterNames{N} info::Vector{N} end {::String,::Int}\n@tree struct SubTreeWithWrongTypeParameterNames{N} info::Vector{N} end {<:AbstractString,<:Number}\nwill get the same error message: \"@tree error: :N and :D are reserved type parameter names.\"\nWhen the formal argument names of keytype and valtype overlap with those of the raw struct type parameters, the duplicates will be considered as the same. For example, the decorated struct SubTreeWithOverlappedParametricFields by the following code\n@tree struct SubTreeWithOverlappedParametricFields{N} info::Vector{N} end {N<:AbstractString,D<:Number}\nonly has two type parameters N<:AbstractString and D<:Number, where the N in the info::Vector{N} is the same N with that in the decorated attribute TREECORE::TreeCore{N,D}. While, if the formal names of keytype and valtype have no intersection with those of the raw struct type parameters, the type parameters of the decorated struct will be just extended by keytype and valtype. For example, the decorated struct SubTreeWithParametricFields by the following code\n@tree struct SubTreeWithParametricFields{T} info::Vector{T} end {N<:AbstractString,D<:Number}\nhave 3 type parameters, T, N<:AbstractString and D<:Number.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Hamiltonian.Utilities.Tree.addnode!-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},N}} where D where N",
    "page": "Tree",
    "title": "Hamiltonian.Utilities.Tree.addnode!",
    "category": "method",
    "text": "addnode!(tree::AbstractTree{N,D},node::N) where {N,D} -> typeof(tree)\naddnode!(tree::AbstractTree{N,D},::Nothing,node::N) where {N,D} -> typeof(tree)\naddnode!(tree::AbstractTree{N,D},parent::N,node::N) where {N,D} -> typeof(tree)\n\nUpdate the structure of a tree by adding a node. When the parent is nothing, the input tree must be empty and the input node becomes the tree\'s root.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Hamiltonian.Utilities.Tree.ancestor-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},N}, Tuple{AbstractTree{N,D},N,Int64}} where D where N",
    "page": "Tree",
    "title": "Hamiltonian.Utilities.Tree.ancestor",
    "category": "method",
    "text": "ancestor(tree::AbstractTree{N,D},node::N,generation::Int=1) where {N,D} -> N\n\nGet the ancestor of a tree\'s node of the n-th generation.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Hamiltonian.Utilities.Tree.children-Tuple{Hamiltonian.Utilities.Tree.AbstractTree}",
    "page": "Tree",
    "title": "Hamiltonian.Utilities.Tree.children",
    "category": "method",
    "text": "children(tree::AbstractTree) -> Vector{keytype(tree)}\nchildren(tree::AbstractTree,::Nothing) -> Vector{keytype(tree)}\nchildren(tree::AbstractTree{N,D},node::N) where {N,D} -> Vector{N}\n\nGet the children of a tree\'s node.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Hamiltonian.Utilities.Tree.deletenode!-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},N}} where D where N",
    "page": "Tree",
    "title": "Hamiltonian.Utilities.Tree.deletenode!",
    "category": "method",
    "text": "deletenode!(tree::AbstractTree{N,D},node::N) where {N,D} -> typeof(tree)\n\nUpdate the structure of a tree by deleting a node.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Hamiltonian.Utilities.Tree.descendants-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},N}, Tuple{AbstractTree{N,D},N,Int64}} where D where N",
    "page": "Tree",
    "title": "Hamiltonian.Utilities.Tree.descendants",
    "category": "method",
    "text": "descendants(tree::AbstractTree{N,D},node::N,generation::Int=1) where {N,D} -> Vector{N}\n\nGet the descendants of a tree\'s node of the n-th generation.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Hamiltonian.Utilities.Tree.empty-Tuple{Hamiltonian.Utilities.Tree.AbstractTree}",
    "page": "Tree",
    "title": "Hamiltonian.Utilities.Tree.empty",
    "category": "method",
    "text": "empty(tree::AbstractTree)\n\nConstruct an empty tree of the same type with the input one.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Hamiltonian.Utilities.Tree.isleaf-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},N}} where D where N",
    "page": "Tree",
    "title": "Hamiltonian.Utilities.Tree.isleaf",
    "category": "method",
    "text": "isleaf(tree::AbstractTree{N,D},node::N) where{N,D} -> Bool\n\nJudge whether a tree\'s node is a leaf (a node without children) or not.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Hamiltonian.Utilities.Tree.leaves-Tuple{Hamiltonian.Utilities.Tree.AbstractTree}",
    "page": "Tree",
    "title": "Hamiltonian.Utilities.Tree.leaves",
    "category": "method",
    "text": "leaves(tree::AbstractTree) -> Vector{keytype(tree)}\n\nGet a tree\'s leaves.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Hamiltonian.Utilities.Tree.level-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},N}} where D where N",
    "page": "Tree",
    "title": "Hamiltonian.Utilities.Tree.level",
    "category": "method",
    "text": "level(tree::AbstractTree{N,D},node::N) where {N,D} -> Int\n\nGet the level of tree\'s node.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Hamiltonian.Utilities.Tree.move!-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},N,N}} where D where N",
    "page": "Tree",
    "title": "Hamiltonian.Utilities.Tree.move!",
    "category": "method",
    "text": "move!(tree::AbstractTree{N,D},node::N,parent::N) where {N,D} -> typeof(tree)\n\nMove a subtree to a new position.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Hamiltonian.Utilities.Tree.parent-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},N}, Tuple{AbstractTree{N,D},N,Union{Nothing, N}}} where D where N",
    "page": "Tree",
    "title": "Hamiltonian.Utilities.Tree.parent",
    "category": "method",
    "text": "parent(tree::AbstractTree{N,D},node::N,superparent::Union{N,Nothing}=nothing) where {N,D} -> Union{N,Nothing}\n\nGet the parent of a tree\'s node. When node is the tree\'s root, return superparent.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Hamiltonian.Utilities.Tree.root-Tuple{Hamiltonian.Utilities.Tree.AbstractTree}",
    "page": "Tree",
    "title": "Hamiltonian.Utilities.Tree.root",
    "category": "method",
    "text": "root(tree::AbstractTree) -> Union{keytype(tree),Nothing}\n\nGet a tree\'s root node.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Hamiltonian.Utilities.Tree.siblings-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},N}} where D where N",
    "page": "Tree",
    "title": "Hamiltonian.Utilities.Tree.siblings",
    "category": "method",
    "text": "siblings(tree::AbstractTree{N,D},node::N) where{N,D} -> Vector{N}\n\nGet the siblings (other nodes sharing the same parent) of a tree\'s node.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Hamiltonian.Utilities.Tree.subtree-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},N}} where D where N",
    "page": "Tree",
    "title": "Hamiltonian.Utilities.Tree.subtree",
    "category": "method",
    "text": "subtree(tree::AbstractTree{N,D},node::N) where{N,D} -> typeof(tree)\n\nGet a subtree whose root is node.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Base.:==-Union{Tuple{TC}, Tuple{TC,TC}} where TC<:Hamiltonian.Utilities.Tree.TreeCore",
    "page": "Tree",
    "title": "Base.:==",
    "category": "method",
    "text": "==(tc1::TC,tc2::TC) where TC<:TreeCore -> Bool\n\nOverloaded == operator.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Base.:==-Union{Tuple{T}, Tuple{T,T}} where T<:Hamiltonian.Utilities.Tree.AbstractTree",
    "page": "Tree",
    "title": "Base.:==",
    "category": "method",
    "text": "==(t1::T,t2::T) where T<:AbstractTree -> Bool\n\nOverloaded == operator.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Base.append!-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},AbstractTree{N,D}}} where D where N",
    "page": "Tree",
    "title": "Base.append!",
    "category": "method",
    "text": "append!(tree::AbstractTree{N,D},subtree::AbstractTree{N,D}) where {N,D} -> typeof(tree)\nappend!(tree::AbstractTree{N,D},node::Union{N,Nothing},subtree::AbstractTree{N,D}) where {N,D} -> typeof(tree)\n\nAppend a subtree to a tree.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Base.delete!-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},N}} where D where N",
    "page": "Tree",
    "title": "Base.delete!",
    "category": "method",
    "text": "delete!(tree::AbstractTree{N,D},node::N) where {N,D} -> typeof(tree)\n\nDelete a node and all its descendants from a tree.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Base.eltype-Tuple{Hamiltonian.Utilities.Tree.AbstractTree}",
    "page": "Tree",
    "title": "Base.eltype",
    "category": "method",
    "text": "eltype(tree::AbstractTree)\neltype(::Type{<:AbstractTree{N,D}}) where {N,D}\n\nGet the eltype of a tree.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Base.empty!-Tuple{Hamiltonian.Utilities.Tree.AbstractTree}",
    "page": "Tree",
    "title": "Base.empty!",
    "category": "method",
    "text": "empty!(tree::AbstractTree) -> typeof(tree)\n\nEmpty a tree.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Base.getindex-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},N}} where D where N",
    "page": "Tree",
    "title": "Base.getindex",
    "category": "method",
    "text": "getindex(tree::AbstractTree{N,D},node::N) where {N,D} -> N\n\nGet the data of a tree\'s node.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Base.haskey-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},N}} where D where N",
    "page": "Tree",
    "title": "Base.haskey",
    "category": "method",
    "text": "haskey(tree::AbstractTree{N,D},node::N) where {N,D} -> Bool\n\nCheck whether a node is in a tree.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Base.keys-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},Val{0}}, Tuple{AbstractTree{N,D},Val{0},Union{Nothing, N}}} where D where N",
    "page": "Tree",
    "title": "Base.keys",
    "category": "method",
    "text": "keys(tree::AbstractTree{N,D},::typeof(treedepth),node::Union{N,Nothing}=tree|>root) where {N,D}\nkeys(tree::AbstractTree{N,D},::typeof(treewidth),node::Union{N,Nothing}=tree|>root) where {N,D}\n\nIterate over a tree\'s nodes starting from a certain node by depth first search or width first search.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Base.keytype-Tuple{Hamiltonian.Utilities.Tree.AbstractTree}",
    "page": "Tree",
    "title": "Base.keytype",
    "category": "method",
    "text": "keytype(tree::AbstractTree)\nkeytype(::Type{<:AbstractTree{N,D}}) where {N,D}\n\nGet a tree\'s node type.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Base.length-Tuple{Hamiltonian.Utilities.Tree.AbstractTree}",
    "page": "Tree",
    "title": "Base.length",
    "category": "method",
    "text": "length(tree::AbstractTree) -> Int\n\nGet the number of a tree\'s nodes.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Base.pairs-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree,Val{0}}, Tuple{AbstractTree,Val{0},Union{Nothing, N}}} where D where N",
    "page": "Tree",
    "title": "Base.pairs",
    "category": "method",
    "text": "pairs(tree::AbstractTree,::typeof(treedepth),node::Union{N,Nothing}=tree|>root) where {N,D}\npairs(tree::AbstractTree,::typeof(treewidth),node::Union{N,Nothing}=tree|>root) where {N,D}\n\nIterate over a tree\'s (node,data) pairs starting from a certain node by depth first search or width first search.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Base.push!-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},N,D}} where D where N",
    "page": "Tree",
    "title": "Base.push!",
    "category": "method",
    "text": "push!(tree::AbstractTree{N,D},node::N,data::D) where {N,D} -> typeof(tree)\npush!(tree::AbstractTree{N,D},parent::Union{N,Nothing},node::N,data::D) where {N,D} -> typeof(tree)\n\nPush a new node to a tree. When parent is nothing, this function set the root node of an empty tree.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Base.setindex!-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},D,N}} where D where N",
    "page": "Tree",
    "title": "Base.setindex!",
    "category": "method",
    "text": "setindex!(tree::AbstractTree{N,D},data::D,node::N) where {N,D}\n\nSet the data of a tree\'s node.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Base.valtype-Tuple{Hamiltonian.Utilities.Tree.AbstractTree}",
    "page": "Tree",
    "title": "Base.valtype",
    "category": "method",
    "text": "valtype(tree::AbstractTree)\nvaltype(::Type{<:AbstractTree{N,D}}) where {N,D}\n\nGet a tree\'s data type.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Base.values-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree,Val{0}}, Tuple{AbstractTree,Val{0},Union{Nothing, N}}} where D where N",
    "page": "Tree",
    "title": "Base.values",
    "category": "method",
    "text": "values(tree::AbstractTree,::typeof(treedepth),node::Union{N,Nothing}=tree|>root) where {N,D}\nvalues(tree::AbstractTree,::typeof(treewidth),node::Union{N,Nothing}=tree|>root) where {N,D}\n\nIterate over a tree\'s data starting from a certain node by depth first search or width first search.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Manual-1",
    "page": "Tree",
    "title": "Manual",
    "category": "section",
    "text": "Modules=[Tree]\nOrder=  [:module,:constant,:type,:macro,:function]"
},

{
    "location": "man/Utilities/NamedVector.html#",
    "page": "Named vector",
    "title": "Named vector",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Utilities.NamedVectorpush!(LOAD_PATH,\"../../../../src/\")\nusing Hamiltonian.Utilities.NamedVector"
},

{
    "location": "man/Utilities/NamedVector.html#Named-vector-1",
    "page": "Named vector",
    "title": "Named vector",
    "category": "section",
    "text": "A named vector is similiar to a named tuple, which associate each of its values with a name. Although the names of a named vector cannot be changed, the values can be modified if needed. In contrast to the predefined NamedTuple in Julia, which employs the names as type parameters, we just implement a named vector as a composite struct equipped with the getindex and setindex! functions, with the fieldnames being its names. This simple implementation makes it possible to define your own concrete named vector with any of your preferred type names, and ensures that all instances of a certain concrete named vector share the same names. Therefore, if you are familiar with Python, you will find that our named vector is more qualified to be the counterpart of the namedtuple in Python than the default Julia implementation. Last but not least important, it is also worth noted that a named vector is not a vector, as is similar to that a named tuple is not a tuple in Julia. This results from our basic expectation that a named vector should be more like a tuple other than a vector so that not all operations valid to vectors are also valid to named vectors."
},

{
    "location": "man/Utilities/NamedVector.html#AbstractNamedVector-1",
    "page": "Named vector",
    "title": "AbstractNamedVector",
    "category": "section",
    "text": "AbstractNamedVector defines the abstract type for all concrete named vectors.Main features include:Values can be accessed or modified either by the . operator or by the [] operator.\nComparisons, such as ≡, ≢, ==, ≠, >, <, ≥, ≤ are supported. Therefore a vector of named vectors can be sorted by the default sort function.\nHash is supported by hash. Therefore, a named vector can be used as the key of a dict or set.\nIteration over its fieldnames is supported by keys, over its values is supported by values, over its field-value pairs is supported by pairs. A reverse iteration is also supported.To subtype it, please note:A concrete type can be either mutable or immutable as you need, which is different from tuples.\nThe fields of a concrete type can be of the same type or not. For the former, we denote the named vector as \"homogeneous\" while for the latter as \"inhomogeneous\". For homogeneous ones, we define a sub abstract type, HomoNamedVector for further optimization of the default methods. See HomoNamedVector below.\nIt is recommended to overload the Base.fieldnames function for concrete subtypes to ensure type stability and improve efficiency, which though is not a necessity. A template for such an overloading is\nBase.fieldnames(Type{<:YourNamedVector})=(:fieldname1,:fieldname2,...)\nFor all concrete subtypes, if inner constructors are defined, the one which has the same interface with the default one must be implemented. Otherwise, some functionalities will not work.\nArithmetic operations, such as +, -, *, /, %, ÷, etc. are NOT supported. However, the function map is implemented, which can help users do the overloadings of these operations.We define a macro @namedvector as the type factory to decorate a \"raw\" struct to be a subtype of AbstractNamedVector. Here, \"raw\" means the struct to be decorated has no explicit supertype other than Any, neither inner constructors as well. For example,@namedvector mutable struct InHomoNV\n    scope::String\n    site::Int\nendThis macro encapsulate the overloading of Base.fieldnames, and you have no need to do this by hand any more."
},

{
    "location": "man/Utilities/NamedVector.html#HomoNamedVector-1",
    "page": "Named vector",
    "title": "HomoNamedVector",
    "category": "section",
    "text": "HomoNamedVector is the subtype of [AbstractNamedVector] that of all its fields share the same type. Compared to AbstractNamedVector, one more default method is implemented with HomoNamedVector, i.e. eltype, which returns the type of its fields. This function ensures the type stability of all the methods that involves an iteration of the field values of a named vector. Therefore, homogeneous named vector are usually more efficient than inhomogeneous ones. Use homogeneous ones as much as possible unless the code efficiency does not matter.To subtype [HomoNamedVector], all the suggestions mentioned in the previous subsection for AbstractNamedVector also applies. A recommended template for a subtype is[mutable] struct YourNamedVector{T} <: HomoNamedVector{T}\n    filed1::T\n    filed2::T\n    ...\nendWe also provide a macro @homonamedvector to help the definition of concrete homogeneous named vector, where you only need specify the type name, field names, data type and optionally whether the subtype is mutable. For example,@homonamedvector HomoNVWithoutParameter (:scope,:site) Int mutable=true\n@homonamedvector HomoNVWithParameter (:scope,:site) (<:Real) mutable=trueThis macro also integrates the Base.fieldnames function, thus its overloading by hand is on longer needed."
},

{
    "location": "man/Utilities/NamedVector.html#Hamiltonian.Utilities.NamedVector.AbstractNamedVector",
    "page": "Named vector",
    "title": "Hamiltonian.Utilities.NamedVector.AbstractNamedVector",
    "category": "type",
    "text": "AbstractNamedVector\n\nAbstract type for all named vectors.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Hamiltonian.Utilities.NamedVector.HomoNamedVector",
    "page": "Named vector",
    "title": "Hamiltonian.Utilities.NamedVector.HomoNamedVector",
    "category": "type",
    "text": "HomoNamedVector{T}\n\nAbstract type for all homogeneous named vectors.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Hamiltonian.Utilities.NamedVector.@homonamedvector",
    "page": "Named vector",
    "title": "Hamiltonian.Utilities.NamedVector.@homonamedvector",
    "category": "macro",
    "text": "@homonamedvector typename fieldnames dtype::Union{Expr,Symbol}=:nothing mutable::Union{Expr,Bool}=false\n\nConstruct a concrete homogeneous named vector with the type name being typename and the fieldnames specified by fieldnames, and optionally, the type parameters specified by dtype.mutable can be used as a keyword argument to determine whether the concrete type is mutable.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Hamiltonian.Utilities.NamedVector.@namedvector-Tuple{Expr}",
    "page": "Named vector",
    "title": "Hamiltonian.Utilities.NamedVector.@namedvector",
    "category": "macro",
    "text": "@namedvector structdef::Expr\n\nDecorate a \"raw\" struct to be a subtype of AbstractNamedVector. Here, \"raw\" means that the input struct has no explicit supertype and no inner constructors.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.:<-Union{Tuple{NV}, Tuple{NV,NV}} where NV<:Hamiltonian.Utilities.NamedVector.AbstractNamedVector",
    "page": "Named vector",
    "title": "Base.:<",
    "category": "method",
    "text": "<(nv1:NV,nv2:NV) where NV<:AbstractNamedVector -> Bool\n\nOverloaded < operator.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.:==-Tuple{Hamiltonian.Utilities.NamedVector.AbstractNamedVector,Hamiltonian.Utilities.NamedVector.AbstractNamedVector}",
    "page": "Named vector",
    "title": "Base.:==",
    "category": "method",
    "text": "==(nv1::AbstractNamedVector,nv2::AbstractNamedVector) -> Bool\n\nOverloaded == operator. Two named vector are equal to each other if and only if their keys as well as their values are equal to each other.\n\nnote: Note\nIt is not necessary for two named vectors to be of the same concrete type to be equal to each other.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.convert-Tuple{Type{Tuple},Hamiltonian.Utilities.NamedVector.AbstractNamedVector}",
    "page": "Named vector",
    "title": "Base.convert",
    "category": "method",
    "text": "convert(::Type{Tuple},nv::AbstractNamedVector) -> Tuple\nconvert(::Type{NV},nv::Tuple) where NV<:AbstractNamedVector -> NV\n\nConvert a named vector to tuple and vice versa.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.eltype-Union{Tuple{Type{#s68} where #s68<:HomoNamedVector{T}}, Tuple{T}} where T",
    "page": "Named vector",
    "title": "Base.eltype",
    "category": "method",
    "text": "eltype(::Type{NV}) where NV<:HomoNamedVector{T} where T\neltype(nv::HomoNamedVector)\n\nGet the type parameter of a concrete HomoNamedVector.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.getindex-Tuple{Hamiltonian.Utilities.NamedVector.AbstractNamedVector,Int64}",
    "page": "Named vector",
    "title": "Base.getindex",
    "category": "method",
    "text": "getindex(nv::AbstractNamedVector,index::Int)\n\nGet the value by the [] syntax.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.hash-Tuple{Hamiltonian.Utilities.NamedVector.AbstractNamedVector,UInt64}",
    "page": "Named vector",
    "title": "Base.hash",
    "category": "method",
    "text": "hash(nv::AbstractNamedVector,h::UInt)\n\nHash a concrete AbstractNamedVector.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.isless-Union{Tuple{NV}, Tuple{NV,NV}} where NV<:Hamiltonian.Utilities.NamedVector.AbstractNamedVector",
    "page": "Named vector",
    "title": "Base.isless",
    "category": "method",
    "text": "isless(nv1::NV,nv2::NV) where NV<:AbstractNamedVector -> Bool\n\nOverloaded isless function.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.iterate",
    "page": "Named vector",
    "title": "Base.iterate",
    "category": "function",
    "text": "iterate(nv::AbstractNamedVector,state=1)\niterate(rv::Iterators.Reverse{<:AbstractNamedVector},state=length(rv.itr))\n\nIterate or reversely iterate over the values of a concrete AbstractNamedVector.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.keys-Tuple{Hamiltonian.Utilities.NamedVector.AbstractNamedVector}",
    "page": "Named vector",
    "title": "Base.keys",
    "category": "method",
    "text": "keys(nv::AbstractNamedVector) -> NTuple(nv|>length,Symbol)\n\nIterate over the names.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.length-Union{Tuple{Type{NV}}, Tuple{NV}} where NV<:Hamiltonian.Utilities.NamedVector.AbstractNamedVector",
    "page": "Named vector",
    "title": "Base.length",
    "category": "method",
    "text": "length(::Type{NV}) where NV<:AbstractNamedVector -> Int\nlength(nv::AbstractNamedVector) -> Int\n\nGet the length of a concrete AbstractNamedVector.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.map-Union{Tuple{NV}, Tuple{Any,Vararg{NV,N} where N}} where NV<:Hamiltonian.Utilities.NamedVector.AbstractNamedVector",
    "page": "Named vector",
    "title": "Base.map",
    "category": "method",
    "text": "map(f,nvs::NV...) where NV<:AbstractNamedVector -> NV\n\nApply function f elementwise on the input named vectors.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.pairs-Tuple{Hamiltonian.Utilities.NamedVector.AbstractNamedVector}",
    "page": "Named vector",
    "title": "Base.pairs",
    "category": "method",
    "text": "pairs(nv::AbstractNamedVector)\n\nIterate over the name=>value pairs.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.replace-Tuple{Hamiltonian.Utilities.NamedVector.AbstractNamedVector}",
    "page": "Named vector",
    "title": "Base.replace",
    "category": "method",
    "text": "replace(nv::AbstractNamedVector;kwargs...) -> typeof(nv)\n\nReturn a copy of a concrete AbstractNamedVector with some of the field values replaced by the keyword arguments.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.setindex!-Tuple{Hamiltonian.Utilities.NamedVector.AbstractNamedVector,Any,Int64}",
    "page": "Named vector",
    "title": "Base.setindex!",
    "category": "method",
    "text": "setindex!(nv::AbstractNamedVector,value,index::Int)\n\nSet the value by the [] syntax if mutable.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.show-Tuple{IO,Hamiltonian.Utilities.NamedVector.AbstractNamedVector}",
    "page": "Named vector",
    "title": "Base.show",
    "category": "method",
    "text": "show(io::IO,nv::AbstractNamedVector)\n\nShow a concrete AbstractNamedVector.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.values-Tuple{Hamiltonian.Utilities.NamedVector.AbstractNamedVector}",
    "page": "Named vector",
    "title": "Base.values",
    "category": "method",
    "text": "values(nv::AbstractNamedVector) -> Tuple\n\nIterate over the values.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.zero-Union{Tuple{Type{NV}}, Tuple{NV}} where NV<:Hamiltonian.Utilities.NamedVector.AbstractNamedVector",
    "page": "Named vector",
    "title": "Base.zero",
    "category": "method",
    "text": "zero(::Type{NV}) where NV<:AbstractNamedVector -> NV\nzero(nv::AbstractNamedVector) -> typeof(nv)\n\nGet a concrete AbstractNamedVector with all values being zero.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Manual-1",
    "page": "Named vector",
    "title": "Manual",
    "category": "section",
    "text": "Modules=[NamedVector]\nOrder=  [:module,:constant,:type,:macro,:function]"
},

{
    "location": "man/Utilities/QuantumNumber.html#",
    "page": "Quantum numbers",
    "title": "Quantum numbers",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Utilities.QuantumNumber"
},

{
    "location": "man/Utilities/QuantumNumber.html#Quantum-numbers-1",
    "page": "Quantum numbers",
    "title": "Quantum numbers",
    "category": "section",
    "text": "Qunatum numbers can be considered as the conserved labels for the bases of a Hilbert space when a quantum system hosts some symmetries. Here we only implement Abelian quantum numbers because non-Abelian ones are far more complicated yet much less used. In practice, quantum numbers can be integers or half integers, therefore, we use real numbers to denote them in this module for simplicity. Independent quantum numbers, such as the particle number and the spin z-component, can coexist at the same time. We use type AbstractQuantumNumber to represent the complete set of independent ones for a single basis of a Hilbert space, and type QuantumNumbers to represent the whole quantum numbers for the total bases."
},

{
    "location": "man/Utilities/QuantumNumber.html#AbstractQuantumNumber-1",
    "page": "Quantum numbers",
    "title": "AbstractQuantumNumber",
    "category": "section",
    "text": "The abstract type for the complete set of independent quantum numbers for a single basis.Main features include:function fieldnames: get the names of the quantum numbers\nfunction periods: get the periods of the quantum numbers\narithmetic operations: +, -, *, ^, ⊕, ⊗\nhashable: concrete instances can be used as keys for a dict or a set\niterable: concrete instances are iterable over their values\ncomparable: two concrete instances can be comparedIn particular, AbstractQuantumNumber <: AbstractNamedVector{Float64}, all features supported by AbstractNamedVector are also available for AbstractQuantumNumber. See also AbstractNamedVector.For convenience, 4 kinds of quantum numbers are predefined in this module, i.e.SQN: for spin z-component reserved systems\nPQN: for particle number reserved systems\nSPQN: for both particle number and spin-z component reserved systems\nZ2QN: for systems with a Z_2 conservation quantum numberUsers who want to define their own Z_N-like quantum numbers must handle the periodicities in the construction function, otherwise, wrong results will be get when arithmetic operations, such as + or -, are involved. It is highly recommended to use the macro @quantumnumber to define your own concrete AbstractQuantumNumbers."
},

{
    "location": "man/Utilities/QuantumNumber.html#QuantumNumbers-1",
    "page": "Quantum numbers",
    "title": "QuantumNumbers",
    "category": "section",
    "text": "The whole quantum numbers for the total bases.By design, a QuantumNumbers{QN} has one type parameter:QN<:AbstractQuantumNumber: the type of the quantum numbers contained in itAnd 3 attributes:form::Char: Its form, whose value must be one of the followings\n\'G\': the general form, which has no restriction for its contents\n\'U\': the unitary form, which requires no duplicates in its contents\n\'C\': the canonical form, which requires not only no duplicates but also accending-order storage in its contents\nUsually, \'G\'-formed and \'U\'-formed QuantumNumberses can be transformed to the corresponding \'C\'-formed ones by the sort function.\ncontents::Vector{QN}: The quantum numbers contained in it. To achieve high efficiency, it is required to be an homogenous array of a certain kind of concrete AbstractQuantumNumber.\nindptr::Vector{Int}: The indptr of the quantum numbers contained in it, which is similar to the colptr attribute of a CSC sparse matrix and records the compression info of its contents.Main features include:function eltype: get the concrete type of the quantum numbers it contains\nindex access: get the contents directly by the getindex function\narithmetic operations: +, -, *, ^, ⊗, ⊕\niterable: various iteration supports, including functions such as iterate, keys, values and pairs\n...For a complete summation of its features, please refer to the manual.For convenience, 5 functions are predefined to generate the QuantumNumbers of common physical systems, i.e.SQNS: a signle spin\nPQNS: a single-particle state with at most N identical particles\nSzPQNS: a single-paritcle state with at most one particle whose spin-z component is Sz\nSPQNS: a single site with internal degrees of freedom that can be ascribed to a spin\nZ2QNS: any Z_2 Hilbert space"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.qnsbruteforce",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.qnsbruteforce",
    "category": "constant",
    "text": "Choice associated with QuantumNumbers, meaning \'by brute force\'.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.qnscontents",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.qnscontents",
    "category": "constant",
    "text": "Choice associated with QuantumNumbers, meaning \'for contents\'.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.qnscounts",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.qnscounts",
    "category": "constant",
    "text": "Choice associated with QuantumNumbers, meaning \'by counts\'.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.qnsexpansion",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.qnsexpansion",
    "category": "constant",
    "text": "Choice associated with QuantumNumbers, meaning \'for expansion\'.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.qnsindices",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.qnsindices",
    "category": "constant",
    "text": "Choice associated with QuantumNumbers, meaning \'for indices\'.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.qnsindptr",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.qnsindptr",
    "category": "constant",
    "text": "Choice associated with QuantumNumbers, meaning \'by indptr\'.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.qnsmontecarlo",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.qnsmontecarlo",
    "category": "constant",
    "text": "Choice associated with quantumnumbers, meaning \'by Monte Carlo\'.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.AbstractQuantumNumber",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.AbstractQuantumNumber",
    "category": "type",
    "text": "Abstract type for all concrete quantum numbers for a single basis.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.PQN",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.PQN",
    "category": "type",
    "text": "PQN(N::Real)\n\nThe concrete AbstractQuantumNumber of a quantum system with particle number N conserved.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.QuantumNumbers",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.QuantumNumbers",
    "category": "type",
    "text": "QuantumNumbers(form::Char,contents::Vector{<:AbstractQuantumNumber},counts::Vector{Int},::typeof(qnscounts))\nQuantumNumbers(form::Char,contents::Vector{<:AbstractQuantumNumber},indptr::Vector{Int},::typeof(qnsindptr))\n\nThe whole quantum numbers of the total bases of a Hilbert space. The default constructors construct a QuantumNumbers from a vector of concrete quantum numbers and an vector containing their counts or indptr.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.QuantumNumbers",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.QuantumNumbers",
    "category": "type",
    "text": "QuantumNumbers(qn::AbstractQuantumNumber,count::Int=1)\n\nConstruct a QuantumNumbers with one unique quantum number which occurs count times.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.QuantumNumbers-Tuple{OrderedCollections.OrderedDict{#s15,Int64} where #s15<:Hamiltonian.Utilities.QuantumNumber.AbstractQuantumNumber}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.QuantumNumbers",
    "category": "method",
    "text": "QuantumNumbers(od::OrderedDict{<:AbstractQuantumNumber,Int})\n\nConstruct a QuantumNumbers from an ordered dict containing concrete quantum numbers and their counts.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.QuantumNumbers-Tuple{OrderedCollections.OrderedDict{#s15,UnitRange{Int64}} where #s15<:Hamiltonian.Utilities.QuantumNumber.AbstractQuantumNumber}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.QuantumNumbers",
    "category": "method",
    "text": "QuantumNumbers(od::OrderedDict{<:AbstractQuantumNumber,UnitRange{Int}})\n\nConstruct a QuantumNumbers from an ordered dict containing concrete quantum numbers and their slices.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.SPQN",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.SPQN",
    "category": "type",
    "text": "SPQN(N::Real,Sz::Real)\n\nThe concrete AbstractQuantumNumber of a quantum system with both particle number N and spin z-component Sz conserved.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.SQN",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.SQN",
    "category": "type",
    "text": "SQN(Sz::Real)\n\nThe concrete AbstractQuantumNumber of a quantum system with spin z-component Sz conserved.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.Z2QN",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.Z2QN",
    "category": "type",
    "text": "Z2QN(N::Real)\n\nThe concrete AbstractQuantumNumber of a quantum system with a Z₂-like conserved quantity.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.@quantumnumber-Tuple{Any,Any,Any}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.@quantumnumber",
    "category": "macro",
    "text": "@quantumnumber typename fieldnames fieldperiods\n\nConstruct a concrete AbstractQuantumNumber with the type name being typename, fieldnames specified by fieldnames and periods specified by fieldperiods.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.:⊕-Union{Tuple{Tuple{Vararg{#s15,N}} where #s15<:AbstractQuantumNumber}, Tuple{N}, Tuple{Tuple{Vararg{#s14,N}} where #s14<:AbstractQuantumNumber,Tuple{Vararg{Int64,N}}}} where N",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.:⊕",
    "category": "method",
    "text": "⊕(qns::NTuple{N,<:AbstractQuantumNumber},signs::NTuple{N,Int}=ntuple(i->1,N)) where N -> QuantumNumbers\n⊕(qnses::NTuple{N,QuantumNumbers{QN}},signs::NTuple{N,Int}=ntuple(i->1,N)) where {N,QN<:AbstractQuantumNumber} -> QuantumNumbers{QN}\n\nGet the direct sum of some AbstractQuantumNumbers or QuantumNumberses.\n\nnote: Note\nPhysically, the direct sum of a couple of AbstractQuantumNumbers or QuantumNumberses is defined by the direct sum of the bases of the Hilbert spaces they represent. Therefore, the input AbstractQuantumNumbers or QuantumNumberses must be homogenous. Inhomogenous \'AbstractQuantumNumber\'s must be direct producted first to ensure homogenity before the direct sum.\nApparently, the dimension of the result equals the summation of those of the inputs, which means, even for AbstractQuantumNumbers, the result will be naturally a QuantumNumbers because the dimension of the result is largeer than 1.\nSigns of AbstractQuantumNumbers or QuantumNumberses can be provided when getting their direct sums.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.:⊗-Union{Tuple{QN}, Tuple{Type{QN},AbstractQuantumNumber,AbstractQuantumNumber}} where QN<:Hamiltonian.Utilities.QuantumNumber.AbstractQuantumNumber",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.:⊗",
    "category": "method",
    "text": "⊗(::Type{QN},qn1::AbstractQuantumNumber,qn2::AbstractQuantumNumber) where QN<:AbstractQuantumNumber -> QN\n⊗(qns::NTuple{N,<:AbstractQuantumNumber},signs::NTuple{N,Int}=ntuple(i->1,N)) where N -> eltype(qns)\n⊗(qnses::NTuple{N,QuantumNumbers{QN}},signs::NTuple{N,Int}=ntuple(i->1,N)) where {N,QN<:AbstractQuantumNumber} -> QuantumNumbers{QN}\n\nGet the direct product of some AbstractQuantumNumbers or QuantumNumberses.\n\nnote: Note\nPhysically, the direct product of a couple of AbstractQuantumNumbers or QuantumNumberses are defined by the direct product of the bases of the Hilbert spaces they represent. Therefore, QuantumNumbers with differenct types or QuantumNumberses with differenct eltypes are allowed to be direct producted in principle. However, for simplicity, we only implement a method which handle the situation of two AbstractQuantumNumbers with differenct types. The type of the result should be provided as the first parameter. Note that in this situation, the fieldnames and periods of the result type must be exactly equal to the flattened fieldnames and periods of the two input AbstractQuantumNumbers, which means, even the order of the input AbstractQuantumNumbers matters.\nApparently, the dimension of the result equals the product of those of the inputs. Therefore, the direct product of AbstractQuantumNumbers is also a AbstractQuantumNumber since its dimension is still one.\nFor other situations except the one mentioned in Note.1, the input AbstractQuantumNumbers or QuantumNumberses must be homogenous. Meanwhile, signs can also be provided for these situations. Note that each quantum number in the contents of the result is obtained by a summation of the corresponding quanum numbers out of the inputs with the correct signs. This is a direct observation of the Abelian nature of our quantum numbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.PQNS-Tuple{Real}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.PQNS",
    "category": "method",
    "text": "PQNS(N::Real) -> QuantumNumbers{PQN}\n\nConstruct the QuantumNumbers of the Hilbert space of a single-particle state with at most N identical particles.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.SPQNS-Tuple{Real}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.SPQNS",
    "category": "method",
    "text": "SPQNS(S::Real) -> QuantumNumbers{SPQN}\n\nConstruct the QuantumNumbers of the Hilbert space of a single site with internal degrees of freedom that can be ascribed to a spin S.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.SQNS-Tuple{Real}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.SQNS",
    "category": "method",
    "text": "SQNS(S::Real) -> QuantumNumbers{SQN}\n\nConstruct the QuantumNumbers of the Hilbert space of a signle spin S.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.SzPQNS-Tuple{Real}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.SzPQNS",
    "category": "method",
    "text": "SzPQNS(Sz::Real) -> QuantumNumbers{SPQN}\n\nConstruct the QuantumNumbers of the Hilbert space of a single-paritcle state with at most one particle whose spin-z component is Sz.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.Z2QNS-Tuple{}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.Z2QNS",
    "category": "method",
    "text": "Z2QNS() -> QuantumNumbers{Z2QN}\n\nConstruct the QuantumNumbers of a Z_2 Hilbert space.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.decompose-Union{Tuple{QN}, Tuple{N}, Tuple{Tuple{Vararg{QuantumNumbers{QN},N}},QN,Tuple{Vararg{Int64,N}},Val{6}}} where QN<:Hamiltonian.Utilities.QuantumNumber.AbstractQuantumNumber where N",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.decompose",
    "category": "method",
    "text": "decompose(qnses::NTuple{N,QuantumNumbers{QN}},target::QN,signs::NTuple{N,Int},::typeof(qnsbruteforce);nmax::Int=20) where {N,QN<:AbstractQuantumNumber} -> Vector{NTuple{N,Int}}\ndecompose(qnses::NTuple{N,QuantumNumbers{QN}},target::QN,signs::NTuple{N,Int},::typeof(qnsmontecarlo);nmax::Int=20) where {N,QN<:AbstractQuantumNumber} -> Vector{NTuple{N,Int}}\n\nFind a couple of decompositions of target with respect to qnses.\n\nnote: Note\nA tuple of integers (i₁,i₂,...) is called a decomposition of a given target with respect to the given qnses if and only if they satisfy the \"decomposition rule\":sum_textj textsignstextjtimestextqnsestextjtexti_textj==texttargetThis equation is in fact a kind of a set of restricted linear Diophantine equations. Indeed, our quantum numbers are always discrete Abelian ones and all instances of a concrete AbstractQuantumNumber forms a module over the ring of integers. Therefore, each quantum number can be represented as a integral multiple of the unit element of the Abelian module, which results in the final reduction of the above equation to a set of linear Diophantine equations. Then finding a decomposition is equivalent to find a solution of the reduced linear Diophantine equations, with the restriction that the quantum numbers constructed from the solution should be in the corresponding qnses. Here we provide two methods to find such decompositions, one is by brute force (qnsbruteforce case), and the other is by Monte Carlo simultatioins (qnsmontecarlo case).\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.expand-Tuple{Hamiltonian.Utilities.QuantumNumber.QuantumNumbers,Val{3}}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.expand",
    "category": "method",
    "text": "expand(qns::QuantumNumbers,::typeof(qnscontents)) -> Vector{qns|>eltype}\nexpand(qns::QuantumNumbers,::typeof(qnsindices)) -> Vector{Int}\n\nExpand the contents (qnscontents case) or indices (qnsindices case) of a QuantumNumbers to the uncompressed form.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.regularize!-Union{Tuple{QN}, Tuple{Type{QN},AbstractArray{Float64,1}}} where QN<:Hamiltonian.Utilities.QuantumNumber.AbstractQuantumNumber",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.regularize!",
    "category": "method",
    "text": "regularize!(::Type{QN},array::AbstractVector{Float}) where QN<:AbstractQuantumNumber -> typeof(array)\nregularize!(::Type{QN},array::AbstractMatrix{Float}) where QN<:AbstractQuantumNumber -> typeof(array)\n\nRegularize the elements of an array in place so that it can represent quantum numbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.regularize-Union{Tuple{QN}, Tuple{Type{QN},Union{AbstractArray{Float64,1}, AbstractArray{Float64,2}}}} where QN<:Hamiltonian.Utilities.QuantumNumber.AbstractQuantumNumber",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.regularize",
    "category": "method",
    "text": "regularize(::Type{QN},array::Union{AbstractVector{Float},AbstractMatrix{Float}}) where {QN<:AbstractQuantumNumber} -> typeof(array)\n\nRegularize the elements of an array and return a copy that can represent quantum numbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.reorder-Tuple{Hamiltonian.Utilities.QuantumNumber.QuantumNumbers,Array{Int64,1},Val{3}}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.reorder",
    "category": "method",
    "text": "reorder(qns::QuantumNumbers,permutation::Vector{Int},::typeof(qnscontents)) -> QuantumNumbers\nreorder(qns::QuantumNumbers,permutation::Vector{Int},::typeof(qnsexpansion)) -> QuantumNumbers\n\nReorder the quantum numbers contained in a QuantumNumbers with a permutation and return the new one. For qnscontents case, the permutation is for the contents of the original QuantumNumbers while for qnsexpansion case, the permutation is for the expansion of the original QuantumNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.subset-Union{Tuple{QN}, Tuple{QuantumNumbers{QN},QN}} where QN<:Hamiltonian.Utilities.QuantumNumber.AbstractQuantumNumber",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.subset",
    "category": "method",
    "text": "subset(qns::QuantumNumbers{QN},target::QN) where QN<:AbstractQuantumNumber -> QuantumNumbers{QN}\nsubset(qns::QuantumNumbers{QN},targets::NTuple{N,QN}) where {N,QN<:AbstractQuantumNumber} -> QuantumNumbers{QN}\n\nFind a subset of a QuantumNumbers by picking out the quantum numbers in targets.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.toordereddict-Tuple{Hamiltonian.Utilities.QuantumNumber.QuantumNumbers,Val{1}}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.toordereddict",
    "category": "method",
    "text": "toordereddict(qns::QuantumNumbers,::typeof(qnsindptr)) -> OrderedDict{qns|>eltype,UnitRange{Int}}\ntoordereddict(qns::QuantumNumbers,::typeof(qnscounts)) -> OrderedDict{qns|>eltype,Int}\n\nConvert a QuantumNumbers to an ordered dict.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.ukron-Union{Tuple{Tuple{Vararg{QuantumNumbers{QN},N}}}, Tuple{QN}, Tuple{N}, Tuple{Tuple{Vararg{QuantumNumbers{QN},N}},Tuple{Vararg{Int64,N}}}} where QN<:Hamiltonian.Utilities.QuantumNumber.AbstractQuantumNumber where N",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.ukron",
    "category": "method",
    "text": "ukron(qnses::NTuple{N,QuantumNumbers{QN}},signs::NTuple{N,Int}=ntuple(i->1,N)) where {N,QN<:AbstractQuantumNumber} -> QuantumNumbers{QN},Dict{QN,Dict{NTuple{N,QN},UnitRange{Int}}}\n\nUnitary Kronecker product of several QuantumNumberses. The product result as well as the records of the product will be returned.\n\nnote: Note\nAll input QuantumNumbers must be \'U\' formed or \'C\' formed.\nSince duplicate quantum number are not allowed in \'U\' formed and \'C\' formed QuantumNumberses, in general, there exists a merge process of duplicate quantum numbers in the product result. Therefore, records are needed to keep track of this process, which will be returned along with the product result. The records are stored in a Dict{QN,Dict{NTuple{N,QN},UnitRange{Int}}} typed dict, in which, for each unduplicate quantum number qn in the product result, there exist a record Dict((qn₁,qn₂,...)=>start:stop,...) telling what quantum numbers (qn₁,qn₂,...) a mereged duplicate qn comes from and what slice start:stop this merged duplicate corresponds.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.dimension-Tuple{Hamiltonian.Utilities.QuantumNumber.QuantumNumbers}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.dimension",
    "category": "method",
    "text": "dimension(qns::QuantumNumbers) -> Int\n\nThe dimension of the Hilbert space a QuantumNumbers represents.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.dimension-Tuple{Type{#s68} where #s68<:Hamiltonian.Utilities.QuantumNumber.AbstractQuantumNumber}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.dimension",
    "category": "method",
    "text": "dimension(::Type{<:AbstractQuantumNumber}) -> Int\ndimension(::AbstractQuantumNumber) -> Int\n\nThe dimension of the Hilbert space a AbstractQuantumNumber represents. Apparently, this is always 1.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Base.:*-Tuple{Hamiltonian.Utilities.QuantumNumber.AbstractQuantumNumber,Integer}",
    "page": "Quantum numbers",
    "title": "Base.:*",
    "category": "method",
    "text": "*(qn::AbstractQuantumNumber,factor::Integer) -> typeof(qn)\n*(factor::Integer,qn::AbstractQuantumNumber) -> typeof(qn)\n*(qns::QuantumNumbers,factor::Integer) -> QuantumNumbers\n*(factor::Integer,qns::QuantumNumbers) -> QuantumNumbers\n\nOverloaded * operator for the multiplication between an integer and a AbstractQuantumNumber or a QuantumNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Base.:+-Tuple{Hamiltonian.Utilities.QuantumNumber.AbstractQuantumNumber}",
    "page": "Quantum numbers",
    "title": "Base.:+",
    "category": "method",
    "text": "+(qn::AbstractQuantumNumber) -> typeof(qn)\n+(qn::QN,qns::QN...) where QN<:AbstractQuantumNumber -> QN\n+(qns::QuantumNumbers) -> QuantumNumbers\n+(qn::QN,qns::QuantumNumbers{QN}) where QN<:AbstractQuantumNumber -> QuantumNumbers{QN}\n+(qns::QuantumNumbers{QN},qn::QN) where QN<:AbstractQuantumNumber -> QuantumNumbers{QN}\n\nOverloaded + operator for AbstractQuantumNumber and QuantumNumbers.\n\nnote: Note\nThe addition between a QuantumNumbers and a AbstractQuantumNumber is just a global shift of the contents of the QuantumNumbers by the AbstractQuantumNumber, therefore, the result is a QuantumNumbers.\n+ cannot be used between two QuantumNumbers because the result is ambiguous. Instead, use ⊕ for direct sum and ⊗ for direct product.\nTo ensure type stability, two AbstractQuantumNumber can be added together if and only if they are of the same type.\nSimilarly, a AbstractQuantumNumber and a QuantumNumbers can be added together if and only if the former\'s type is the same with the latter\'s eltype.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Base.:--Tuple{Hamiltonian.Utilities.QuantumNumber.AbstractQuantumNumber}",
    "page": "Quantum numbers",
    "title": "Base.:-",
    "category": "method",
    "text": "-(qn::AbstractQuantumNumber) -> typeof(qn)\n-(qn1::QN,qn2::QN) where QN<:AbstractQuantumNumber -> QN\n-(qns::QuantumNumbers) -> QuantumNumbers\n-(qn::QN,qns::QuantumNumbers{QN}) where QN<:AbstractQuantumNumber -> QuantumNumbers{QN}\n-(qns::QuantumNumbers{QN},qn::QN) where QN<:AbstractQuantumNumber -> QuantumNumbers{QN}\n\nOverloaded - operator for AbstractQuantumNumber and QuantumNumbers.\n\nnote: Note\nThe subtraction between a QuantumNumbers and a AbstractQuantumNumber is just a global shift of the contents of the QuantumNumbers by the AbstractQuantumNumber, therefore, the result is a QuantumNumbers.\n- cannot be used between two QuantumNumbers because the result is ambiguous. Instead, use ⊕ with signs for direct sum and ⊗ with signs for direct product.\nTo ensure type stability, a AbstractQuantumNumber can be subtracted by another AbstractQuantumNumber if and only if they are of the same type.\nSimilarly, a AbstractQuantumNumber can be subtracted by a QuantumNumbers or vice versa if and only if the former\'s type is the same with the latter\'s eltype.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Base.:==-Tuple{Hamiltonian.Utilities.QuantumNumber.QuantumNumbers,Hamiltonian.Utilities.QuantumNumber.QuantumNumbers}",
    "page": "Quantum numbers",
    "title": "Base.:==",
    "category": "method",
    "text": "==(qns1::QuantumNumbers,qns2::QuantumNumbers) -> Bool\n\nOverloaded == operator. Two QuantumNumberses are equal to each other if and only if both their contentses and indptrs are elementwise equal to each other.\n\nnote: Note\nIt is not necessary for two QuantumNumberses to have the same eltype nor the same form to be equal to each other.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Base.:^-Tuple{Hamiltonian.Utilities.QuantumNumber.AbstractQuantumNumber,Integer}",
    "page": "Quantum numbers",
    "title": "Base.:^",
    "category": "method",
    "text": "^(qn::AbstractQuantumNumber,factor::Integer) -> typeof(qn)\n^(qns::QuantumNumbers,factor::Integer) -> QuantumNumbers\n\nOverloaded ^ operator for AbstractQuantumNumber and QuantumNumbers. This operation translates into the direct product ⊗ of factor copies of qn or qns.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Base.eltype-Union{Tuple{Type{#s68} where #s68<:QuantumNumbers{QN}}, Tuple{QN}} where QN",
    "page": "Quantum numbers",
    "title": "Base.eltype",
    "category": "method",
    "text": "eltype(::Type{<:QuantumNumbers{QN}}) where QN\neltype(qns::QuantumNumbers)\n\nGet the type of the concrete AbstractQuantumNumber contained in a QuantumNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Base.findall-Union{Tuple{QN}, Tuple{QuantumNumbers{QN},QN,Union{Val{3}, Val{4}}}} where QN<:Hamiltonian.Utilities.QuantumNumber.AbstractQuantumNumber",
    "page": "Quantum numbers",
    "title": "Base.findall",
    "category": "method",
    "text": "findall(qns::QuantumNumbers{QN},target::QN,choice::Union{typeof(qnscontents),typeof(qnsexpansion)}) where QN<:AbstractQuantumNumber -> Vector{Int}\nfindall(qns::QuantumNumbers{QN},targets::NTuple{N,QN},::typeof(qnscontents)) where {N,QN<:AbstractQuantumNumber} -> Vector{Int}\nfindall(qns::QuantumNumbers{QN},targets::NTuple{N,QN},::typeof(qnsexpansion)) where {N,QN<:AbstractQuantumNumber} -> Vector{Int}\n\nFind all the indices of the target quantum numbers in the contents (qnscontents case) or the expansion (qnsexpansion case) of a QuantumNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Base.getindex-Tuple{Hamiltonian.Utilities.QuantumNumber.QuantumNumbers,Int64}",
    "page": "Quantum numbers",
    "title": "Base.getindex",
    "category": "method",
    "text": "getindex(qns::QuantumNumbers,index::Int) -> eltype(qns)\ngetindex(qns::QuantumNumbers,slice::UnitRange{Int}) -> QuantumNumbers\ngetindex(qns::QuantumNumbers,indices::Vector{Int}) -> QuantumNumbers\n\nOverloaded [] operator.\n\nnote: Note\nFor a QuantumNumbers, all these getindex functions act on its contents, i.e. its compressed data, but not on its expansion, i.e. the uncompressed data. This definition is consistent with the length function.\nWhen the index is an integer, the result is a AbstractQuantumNumber, while when the index is a unit range or a vector of intgers, the result is a QuantumNumbers. The logic is quite reasonable because such behaviors are much alike to those of a vector container.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Base.iterate",
    "page": "Quantum numbers",
    "title": "Base.iterate",
    "category": "function",
    "text": "iterate(qns::QuantumNumbers,state::Int=1)\niterate(rv::Iterators.Reverse{<:QuantumNumbers},state::Int=length(rv.itr,false))\n\nIterate or reversely iterate over the concrete AbstractQuantumNumbers contained in a QuantumNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Base.keys-Tuple{Hamiltonian.Utilities.QuantumNumber.QuantumNumbers}",
    "page": "Quantum numbers",
    "title": "Base.keys",
    "category": "method",
    "text": "keys(qns::QuantumNumbers) -> Vector{qns|>eltype}\n\nIterate over the concrete AbstractQuantumNumbers contained in a QuantumNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Base.kron-Union{Tuple{QN}, Tuple{Type{QN},AbstractQuantumNumber,AbstractQuantumNumber}} where QN<:Hamiltonian.Utilities.QuantumNumber.AbstractQuantumNumber",
    "page": "Quantum numbers",
    "title": "Base.kron",
    "category": "method",
    "text": "kron(::Type{QN},qn1::AbstractQuantumNumber,qn2::AbstractQuantumNumber) where QN<:AbstractQuantumNumber -> QN\nkron(qns::NTuple{N,<:AbstractQuantumNumber},signs::NTuple{N,Int}=ntuple(i->1,N)) where N -> eltype(qns)\nkron(qnses::NTuple{N,QuantumNumbers{QN}},signs::NTuple{N,Int}=ntuple(i->1,N)) where {N,QN<:AbstractQuantumNumber} -> QuantumNumbers{QN}\n\nKronecker product of some AbstractQuantumNumbers or QuantumNumberses. This is defined to be equivalent to the direct product ⊗.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Base.length-Tuple{Hamiltonian.Utilities.QuantumNumber.QuantumNumbers}",
    "page": "Quantum numbers",
    "title": "Base.length",
    "category": "method",
    "text": "length(qns::QuantumNumbers) -> Int\n\nGet the number of unduplicate qunatum numbers in the QuantumNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Base.pairs-Tuple{Hamiltonian.Utilities.QuantumNumber.QuantumNumbers,Union{Val{1}, Val{2}}}",
    "page": "Quantum numbers",
    "title": "Base.pairs",
    "category": "method",
    "text": "pairs(qns::QuantumNumbers,choice::Union{typeof(qnsindptr),typeof(qnscounts)})\n\nIterate over the AbstractQuantumNumber=>slice or AbstractQuantumNumber=>count pairs.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Base.show-Tuple{IO,Hamiltonian.Utilities.QuantumNumber.QuantumNumbers}",
    "page": "Quantum numbers",
    "title": "Base.show",
    "category": "method",
    "text": "show(io::IO,qns::QuantumNumbers)\n\nShow a QuantumNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Base.sort-Tuple{Hamiltonian.Utilities.QuantumNumber.QuantumNumbers}",
    "page": "Quantum numbers",
    "title": "Base.sort",
    "category": "method",
    "text": "sort(qns::QuantumNumbers) -> QuantumNumbers,Vector{Int}\n\nSort the quantum numbers of a AbstractQuantumNumber, return the sorted AbstractQuantumNumber and the permutation array that sorts the expansion of the original QuantumNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Base.string-Tuple{Hamiltonian.Utilities.QuantumNumber.QuantumNumbers}",
    "page": "Quantum numbers",
    "title": "Base.string",
    "category": "method",
    "text": "string(qns::QuantumNumbers) -> String\n\nConvert a QuantumNumbers to string.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Base.values-Tuple{Hamiltonian.Utilities.QuantumNumber.QuantumNumbers,Val{1}}",
    "page": "Quantum numbers",
    "title": "Base.values",
    "category": "method",
    "text": "values(qns::QuantumNumbers,::typeof(qnsindptr))\nvalues(qns::QuantumNumbers,::typeof(qnscounts))\n\nIterate over the slices/counts of the QuantumNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#qnmanual-1",
    "page": "Quantum numbers",
    "title": "Manual",
    "category": "section",
    "text": "Modules=[QuantumNumber]\nOrder=  [:module,:constant,:type,:macro,:function]"
},

]}
