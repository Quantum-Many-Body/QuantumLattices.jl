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
    "location": "tutorial/Unitcell Description.html#",
    "page": "Unitcell Description",
    "title": "Unitcell Description",
    "category": "page",
    "text": ""
},

{
    "location": "tutorial/Unitcell Description.html#Unitcell-Description-1",
    "page": "Unitcell Description",
    "title": "Unitcell Description",
    "category": "section",
    "text": "This is the basic idea behind our framework, i.e. we only describe how the quantum lattice system looks like in a unitcell. To achieve this goal, we divide the problem into three steps:STEP 1: define the unitcell of the lattice;\nSTEP 2: define the local internal Hilbert spaces of the system on the unitcell;\nSTEP 3: define the terms connecting local Hilbert spaces on the sites of the unitcell."
},

{
    "location": "tutorial/Engine App Interface.html#",
    "page": "Engine App Interface",
    "title": "Engine App Interface",
    "category": "page",
    "text": ""
},

{
    "location": "tutorial/Engine App Interface.html#Engine-App-Interface-1",
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
    "text": "Utilities of the Hamiltonian package.Pages=  [\n        \"Factory.md\",\n        \"Tree.md\",\n        \"NamedVector.md\",\n        \"GoodQuantumNumber.md\",\n        ]\nDepth=2"
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
    "text": "The aim of Factory is to provide tools to hack into Julia codes without knowing the details of their abstract syntax trees, so that users can somewhat manipulate the existing codes, modify them and generate new ones. In particular, Factory in this module means the representation of certain blocks of Julia codes by a usual Julia struct. This representation is much easier to comprehend than the canonical Expr representation. In general, we propose two basic requirements that any factory must satisfy:Besides the default constructor, a concrete factory can be constructed from any legal Expr expression that represents the block of codes it aims to represent;\nThe canonical Expr expression of the block of codes that a concrete factory represents can be obtained by \"calling\" the factory itself.These two requirements also defines the basic interfaces to interact with factories.Out of practical purposes, we only implemente 6 kinds of factories, i.e. Argument, Parameter, Field, Block, FunctionFactory and TypeFactory, which represent a function argument, a method or type parameter, a struct field, a begin ... end block, a function itself and a struct itself, respectively. Some of the basic methods making the above requirements fulfilled with these types are based on the powerful functions defined in MacroTools."
},

{
    "location": "man/Utilities/Factory.html#Argument-1",
    "page": "Factory",
    "title": "Argument",
    "category": "section",
    "text": "An Argument has 4 attributes:name::Symbol: the name of the argument\ntype::Union{Symbol,Expr}: the type of the argument\nslurp::Bool: whether the argument should be expanded by ...\ndefault::Any: the default value of the argument, nothing for those with no default valuesAll valid expressions representing the arguments of functions can be passed to the constructor:Argument(:arg)\nArgument(:(arg::ArgType))\nArgument(:(arg::ArgType...))\nArgument(:(arg::ArgType=default))Or you can use the macro @argument for construction directly from an argument declaration:@argument arg::ArgType=defaultThe construction from such expressions is based on the the MacroTools.splitarg function. On the other hand, calling an instance of Argument will get the corresponding Expr expression, e.g.,julia> Argument(:arg,:ArgType,false,:default)()\n:(arg::ArgType=default)This feature is based on the MacroTools.combinearg function."
},

{
    "location": "man/Utilities/Factory.html#Parameter-1",
    "page": "Factory",
    "title": "Parameter",
    "category": "section",
    "text": "A Parameter has 2 attributes:name::Symbol: the name of the parameter\ntype::Union{Symbol,Expr}: the type of the parameterAll expressions that represent a type parameter or a method parameter are allowed to be passed to the constructor:Parameter(:Int)\nParameter(:T)\nParameter(:(<:Number))\nParameter(:(T<:Number))The macro @parameter completes the construction directly from a parameter declaration:@parameter T<:Number"
},

{
    "location": "man/Utilities/Factory.html#Field-1",
    "page": "Factory",
    "title": "Field",
    "category": "section",
    "text": "A Field has 2 attributes:name::Symbol: the name of the field\ntype::Union{Symbol,Expr}: the type of the fieldLegal expressions can be used to construct a Field instance by its constructor:Field(:field)\nField(:(field::FieldType))\nField(:(field::ParametricType{Parameter}))The macro @field is also provided to help the construction directly from a field declaration:@field field::FieldTypeThe construction from these expressions is based on the MacroTools.splitarg function and the convertion to these expressions is based on the MacroTools.combinefield function."
},

{
    "location": "man/Utilities/Factory.html#Block-1",
    "page": "Factory",
    "title": "Block",
    "category": "section",
    "text": "A Block has only one attribute:body::Vector{Any}: the body of the begin ... end blockAny expression can be passed to the constructor of Block:Block(:(x=1))\nBlock(:(x=1;y=2))\nBlock(:(begin x=1 end))\nBlock(quote\n        x=1\n        y=2\n    end)Or you can construct a Block instance directly from any code by the macro @block:@block x=1 y=2The body of a block can also be extended by the push! function or the @push! macro.note: Note\nThe body of a Block is somewhat \"flattened\", i.e. it contains no begin ... end blocks. During the initialization, any such input block will be unblocked and added to the body part by part. So is the push! and @push! processes.\nAll LineNumberNodes generated by the input codes will also be included in the block\'s body. However, you can use rmlines! or @rmlines! to remove them from the body of an existing Block, or use rmlines or @rmlines to get a copy with them removed in the body."
},

{
    "location": "man/Utilities/Factory.html#FunctionFactory-1",
    "page": "Factory",
    "title": "FunctionFactory",
    "category": "section",
    "text": "A FunctionFactory has 6 attributes:name::Symbol: the name of the function\nargs::Vector{Argument}: the positional arguments of the function\nkwargs::Vector{Argument}: the keyword arguments of the function\nrtype::Any: the return type of the function\nparams::Vector{Parameter}: the method parameters specified by the where keyword\nbody::Expr: the body of the functionAll expressions that represent functions are allowed to be passed to the constructor:FunctionFactory(:(f()=nothing))\nFunctionFactory(:(f(x)=x))\nFunctionFactory(:(f(x::Int,y::Int;choice::Function=sum)=choice(x,y)))\nFunctionFactory(:(f(x::T,y::T;choice::Function=sum) where T<:Number=choice(x,y)))\nFunctionFactory(:((f(x::T,y::T;choice::Function=sum)::T) where T<:Number=choice(x,y)))\nFunctionFactory(:(\n    function (f(x::T,y::T;choice::Function=sum)::T) where T<:Number\n        choice(x,y)\n    end\n))\nFunctionFactory(quote\n    function (f(x::T,y::T;choice::Function=sum)::T) where T<:Number\n        choice(x,y)\n    end\nend)Similarly, an instance can also be constructed from the macro @functionfactory:@functionfactory (f(x::T,y::T;choice::Function=sum)::T) where T<:Number=choice(x,y)The construction from and the convertion to such expressions are based on the MacroTools.splitdef and MacroTools.combinedef functions, respectively.note: Note\nBecause the form f{T}(x::T,y::T;choice::Function=sum) has no longer been supported since Julia 0.7, the entry :params in the returned dict by MacroTools.splitarg is always missing. Therefore, we abandon its corresponding field in FunctionFactory but use the attribute :params to denote the :whereparams entry.Other features include:Positional arguments can be added by addargs! or @addargs!\nKeyword arguments can be added by addkwargs! or @addkwargs!\nMethod parameters can be added by addparams! or @addparams!\nBody can be extended by extendbody! or @extendbody!"
},

{
    "location": "man/Utilities/Factory.html#TypeFactory-1",
    "page": "Factory",
    "title": "TypeFactory",
    "category": "section",
    "text": "A TypeFactory has 6 attributes:name::Symbol: the name of the struct\nmutable::Bool: whether or not the struct is mutable\nparams::Vector{Parameter}: the type parameters of the struct\nsupertype::Any: the supertype of the struct\nfields::Vector{Field}: the fields of the struct\nconstructors::Vector{FunctionFactory}: the inner constructors of the structAny expression representing valid struct definitions can be passed to the constructor:TypeFactory(:(struct StructName end))\nTypeFactory(:(struct StructName{T} end))\nTypeFactory(:(struct Child{T} <: Parent{T} where T end))\nTypeFactory(:(\n    struct Child{T} <: Parent{T} where T<:TypeParameter\n        field1::T\n        field2::T\n    end\n))\nTypeFactory(quote\n    struct Child{T} <: Parent{T} where T<:TypeParameter\n        field1::T\n        field2::T\n    end\nend)Also, the macro @typefactory supports the construction directly from a type definition:@typefactory struct Child{T} <: Parent{T} where T<:TypeParameter\n                field1::T\n                field2::T\n                Child(field1::T,field2::T=zero(T)) where T=new{T}(field1,field2)\n            endThe construction from these expressions is based on the MacroTools.splitstructdef function. Meanwhile, the convertion to the corresponding expression from a TypeFactory is based on the MacroTools.combinestructdef function.Other features include:Fields can be added by addfields! or @addfields!\nType parameters can be added by addparams! or @addparams!\nInner constructors can be added by addconstructors! or @addconstructors!"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.Argument",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.Argument",
    "category": "type",
    "text": "Argument(name::Symbol,type::FExpr,slurp::Bool,default::Any)\nArgument(name::Symbol;type::FExpr=:Any,slurp::Bool=false,default::Any=nothing)\nArgument(expr::Expr)\n\nThe struct to describe a argument of a function.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.Argument-Tuple{}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.Argument",
    "category": "method",
    "text": "(a::Argument)()\n\nConvert an Argument to the Expr representation of the argument it describes.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.Block",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.Block",
    "category": "type",
    "text": "Block(parts::FExpr...)\n\nThe struct to describe a begin ... end block.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.Block-Tuple{}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.Block",
    "category": "method",
    "text": "(b::Block)()\n\nConvert a Block to the Expr representation of the begin ... end block it describes.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.Field",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.Field",
    "category": "type",
    "text": "Field(name::Symbol,type::FExpr)\nField(name::Symbol;type::FExpr=:Any)\nField(expr::Expr)\n\nThe struct to describe a field of a struct.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.Field-Tuple{}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.Field",
    "category": "method",
    "text": "(f::Field)()\n\nConvert a Field to the Expr representation of the field it describes.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.FunctionFactory",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.FunctionFactory",
    "category": "type",
    "text": "FunctionFactory(name::Symbol,args::Vector{Argument},kwargs::Vector{Argument},rtype::FExpr,params::Vector{Parameter},body::Block)\nFunctionFactory(    name::Symbol;\n                    args::Vector{Argument}=Argument[],\n                    kwargs::Vector{Argument}=Argument[],\n                    rtype::FExpr=:Any,\n                    params::Vector{Parameter}=Parameter[],\n                    body::Block=Block()\n                    )\nFunctionFactory(expr::Expr)\n\nThe struct to describe a function.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.FunctionFactory-Tuple{}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.FunctionFactory",
    "category": "method",
    "text": "(ff::FunctionFactory)()\n\nConvert a FunctionFactory to the Expr representation of the function it describes.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.Parameter",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.Parameter",
    "category": "type",
    "text": "Parameter(name::Symbol,type::FExpr)\nParameter(name::Symbol;type::FExpr=:Any)\nParameter(expr::FExpr)\n\nThe struct to describe a parameter of a function or a type.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.Parameter-Tuple{}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.Parameter",
    "category": "method",
    "text": "(p::Parameter)()\n\nConvert a Parameter to the Expr representation of the parameter it describes.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.TypeFactory",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.TypeFactory",
    "category": "type",
    "text": "TypeFactory(name::Symbol,mutable::Bool,params::Vector{Parameter},supertype::FExpr,fields::Vector{Field},constructors::Vector{FunctionFactory})\nTypeFactory(    name::Symbol;\n                mutable::Bool=false,\n                params::Vector{Parameter}=Parameter[],\n                supertype::FExpr=:Any,\n                fields::Vector{Field}=Field[],\n                constructors::Vector{FunctionFactory}=FunctionFactory[]\n                )\nTypeFactory(expr::Expr)\n\nThe struct to describe a struct.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.TypeFactory-Tuple{}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.TypeFactory",
    "category": "method",
    "text": "(tf::TypeFactory)()\n\nConvert a TypeFactory to the Expr representation of the struct it describes.\n\n\n\n\n\n"
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
    "text": "@functionfactory expr::Expr\n\nConstruct a FunctionFactory directly from a function definition.\n\n\n\n\n\n"
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
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.addargs!-Tuple{Hamiltonian.Utilities.Factory.FunctionFactory,Vararg{Union{Expr, Symbol},N} where N}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.addargs!",
    "category": "method",
    "text": "addargs!(ff::FunctionFactory,args::FExpr...)\naddargs!(ff::FunctionFactory,args::Argument...)\n\nAdd a couple of positional arguments to a function factory.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.addconstructors!-Tuple{Hamiltonian.Utilities.Factory.TypeFactory,Vararg{Expr,N} where N}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.addconstructors!",
    "category": "method",
    "text": "addconstructors!(tf::TypeFactory,constructors::Expr...)\naddconstructors!(tf::TypeFactory,constructors::FunctionFactory...)\n\nAdd a couple of constructors to a type factory.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.addfields!-Tuple{Hamiltonian.Utilities.Factory.TypeFactory,Vararg{Union{Expr, Symbol},N} where N}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.addfields!",
    "category": "method",
    "text": "addfields!(tf::TypeFactory,fields::FExpr...)\naddfields!(tf::TypeFactory,fields::Field...)\n\nAdd a couple of fields to a type factory.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.addkwargs!-Tuple{Hamiltonian.Utilities.Factory.FunctionFactory,Vararg{Union{Expr, Symbol},N} where N}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.addkwargs!",
    "category": "method",
    "text": "addkwargs!(ff::FunctionFactory,kwargs::FExpr...)\naddkwargs!(ff::FunctionFactory,kwargs::Argument...)\n\nAdd a couple of keyword arguments to a function factory.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.addparams!-Tuple{Union{FunctionFactory, TypeFactory},Vararg{Union{Expr, Symbol},N} where N}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.addparams!",
    "category": "method",
    "text": "addparams!(f::Union{FunctionFactory,TypeFactory},params::FExpr...)\naddparams!(f::Union{FunctionFactory,TypeFactory},params::Parameter...)\n\nAdd a couple of method parameters to a function factory or a type factory.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.extendbody!-Tuple{Hamiltonian.Utilities.Factory.FunctionFactory,Vararg{Union{Expr, Symbol},N} where N}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.extendbody!",
    "category": "method",
    "text": "extendbody!(ff::FunctionFactory,parts::FExpr...)\nextendbody!(ff::FunctionFactory,parts::Block...)\n\nExtend the body of a function factory.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.rmlines!-Tuple{Hamiltonian.Utilities.Factory.Block}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.rmlines!",
    "category": "method",
    "text": "rmlines!(b::Block)\n\nRemove line number nodes in the body of a block.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#MacroTools.rmlines-Tuple{Hamiltonian.Utilities.Factory.Block}",
    "page": "Factory",
    "title": "MacroTools.rmlines",
    "category": "method",
    "text": "rmlines(b::Block)\n\nReturn a copy of a block with the line number nodes removed.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.FExpr",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.FExpr",
    "category": "constant",
    "text": "Factory expression types.\n\n\n\n\n\n"
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
    "text": "==(f1::F,f2::F) where F<:AbstractFactory\n\nOverloaded == operator.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Base.push!-Tuple{Hamiltonian.Utilities.Factory.Block,Vararg{Union{Expr, Symbol},N} where N}",
    "page": "Factory",
    "title": "Base.push!",
    "category": "method",
    "text": "push!(b::Block,parts::FExpr...)\npush!(b::Block,parts::Block...)\n\nPush other parts into the body of a block.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Base.replace-Tuple{Hamiltonian.Utilities.Factory.AbstractFactory}",
    "page": "Factory",
    "title": "Base.replace",
    "category": "method",
    "text": "replace(f::AbstractFactory;kwargs...)\n\nReturn a copy of a concrete AbstractFactory with some of the field values replaced by the keyword arguments.\n\n\n\n\n\n"
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
    "text": "CurrentModule=Hamiltonian.Utilities.Tree"
},

{
    "location": "man/Utilities/Tree.html#Tree-1",
    "page": "Tree",
    "title": "Tree",
    "category": "section",
    "text": ""
},

{
    "location": "man/Utilities/Tree.html#AbstractTree-1",
    "page": "Tree",
    "title": "AbstractTree",
    "category": "section",
    "text": ""
},

{
    "location": "man/Utilities/Tree.html#Hamiltonian.Utilities.Tree.AbstractTree",
    "page": "Tree",
    "title": "Hamiltonian.Utilities.Tree.AbstractTree",
    "category": "type",
    "text": "AbstractTree{Node,Data}\n\nAbstract type for all concrete trees.\n\n\n\n\n\n"
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
    "text": "CurrentModule=Hamiltonian.Utilities.NamedVector"
},

{
    "location": "man/Utilities/NamedVector.html#Named-vector-1",
    "page": "Named vector",
    "title": "Named vector",
    "category": "section",
    "text": "A named vector is similiar to a named tuple, which associate each of its values with a name. Although the names of a named vector cannot be changed, the values can be modified if needed. In contrast to the predefined NamedTuple in Julia, which employs the names as type parameters, we just implement a named vector as a composite struct equipped with the getindex and setindex! functions, with the fieldnames being its names. This simple implementation makes it possible to define your own concrete named vector with any of your preferred type names, and ensures that all instances of a certain concrete named vector share the same names. Therefore, if you are familiar with Python, you will find that our named vector is more qualified to be the counterpart of the namedtuple in Python than the default Julia implementation. Specifically, we also define a macro @namedvector as the type factory to help users to define their own concrete named vectors. Last but not least important, it is also worth noted that a named vector is not a vector, as is similar to that a named tuple is not a tuple in Julia. This results from our basic expectation that a named vector should be more like a tuple other than a vector so that not all operations valid to vectors are also valid to named vectors."
},

{
    "location": "man/Utilities/NamedVector.html#AbstractNamedVector-1",
    "page": "Named vector",
    "title": "AbstractNamedVector",
    "category": "section",
    "text": "AbstractNamedVector defines the abstract type for all concrete named vectors.Main features include:Values can be accessed or modified either by the . operator or by the [] operator.\nComparisons, such as ≡, ≢, ==, ≠, >, <, ≥, ≤ are supported. Therefore a vector of named vectors can be sorted by the default sort function.\nHash is supported by hash. Therefore, a named vector can be used as the key of a dict or set.\nIteration over its fieldnames is supported by keys, over its values is supported by values, over its field-value pairs is supported by pairs. A reverse iteration is also supported.To subtype it, please note:A concrete type can be either mutable or immutable as you need, but all its fields should be of the same type. A recommended template for the subtype is\n[mutable] struct YourNamedVector{T} <: AbstractNamedVector{T}\n    filedname1::T\n    filedname2::T\n    ...\nend\nIt is recommended to overload the Base.fieldnames function for concrete subtypes to ensure type stability and improve efficiency, which though is not a necessity. A template for such an overloading is\nBase.fieldnames(Type{YourNamedVector})=(:fieldname1,:fieldname2,...)\nFor all concrete subtypes, if inner constructors are defined, the one which has the same interface with the default one must be implemented. Otherwise, some functionalities will not work.\nArithmetic operations, such as +, -, *, /, %, ÷, etc. are NOT supported. However, an efficient map function  is implemented, which can help users do the overloadings of these operations."
},

{
    "location": "man/Utilities/NamedVector.html#Hamiltonian.Utilities.NamedVector.AbstractNamedVector",
    "page": "Named vector",
    "title": "Hamiltonian.Utilities.NamedVector.AbstractNamedVector",
    "category": "type",
    "text": "AbstractNamedVector{T}\n\nAbstract type for all concrete named vectors.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Hamiltonian.Utilities.NamedVector.@namedvector",
    "page": "Named vector",
    "title": "Hamiltonian.Utilities.NamedVector.@namedvector",
    "category": "macro",
    "text": "@namedvector mutableornot::Bool typename fieldnames dtype::Union{Expr,Symbol}=:nothing supertypename=:AbstractNamedVector\n\nConstruct a mutable or immutable concrete named vector with the type name being typename and the fieldnames specified by fieldnames, and optionally, the type parameters specified by dtype and the supertype specified by supertypename.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Core.Type-Union{Tuple{Tuple{Vararg{T,N}}}, Tuple{T}, Tuple{N}, Tuple{NV}} where T where N where NV<:Hamiltonian.Utilities.NamedVector.AbstractNamedVector",
    "page": "Named vector",
    "title": "Core.Type",
    "category": "method",
    "text": "(::Type{NV})(values::NTuple{N,T}) where {NV<:AbstractNamedVector,N,T}\n\nConstruct a concrete named vector by a tuple.\n\n\n\n\n\n"
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
    "text": "convert(::Type{Tuple},nv::AbstractNamedVector) -> NTuple{nv|>length,nv|>eltype}\nconvert(::Type{NTuple},nv::AbstractNamedVector) -> NTuple{nv|>length,nv|>eltype}\nconvert(::Type{NTuple{N,T}},nv::AbstractNamedVector{T}) where {N,T} -> NTuple{nv|>length,nv|>eltype}\n\nConvert a named vector to tuple.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.eltype-Union{Tuple{Type{#s68} where #s68<:AbstractNamedVector{T}}, Tuple{T}} where T",
    "page": "Named vector",
    "title": "Base.eltype",
    "category": "method",
    "text": "eltype(::Type{NV}) where NV<:AbstractNamedVector{T} where T\neltype(nv::AbstractNamedVector)\n\nGet the type parameter of a concrete AbstractNamedVector.\n\n\n\n\n\n"
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
    "text": "values(nv::AbstractNamedVector) -> NTuple{nv|>length,nv|>eltype}\n\nIterate over the values.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.zero-Union{Tuple{Type{NV}}, Tuple{NV}} where NV<:Hamiltonian.Utilities.NamedVector.AbstractNamedVector",
    "page": "Named vector",
    "title": "Base.zero",
    "category": "method",
    "text": "zero(::Type{NV}) where NV<:AbstractNamedVector\nzero(nv::AbstractNamedVector)\n\nGet a concrete AbstractNamedVector with all values being zero.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Manual-1",
    "page": "Named vector",
    "title": "Manual",
    "category": "section",
    "text": "Modules=[NamedVector]\nOrder=  [:module,:constant,:type,:macro,:function]"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#",
    "page": "Good quantum numbers",
    "title": "Good quantum numbers",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Utilities.GoodQuantumNumber"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Good-quantum-numbers-1",
    "page": "Good quantum numbers",
    "title": "Good quantum numbers",
    "category": "section",
    "text": "Good qunatum numbers can be considered as the conserved labels for the bases of a Hilbert space when a quantum system hosts some symmetries. Here we only implement Abelian good quantum numbers because non-Abelian ones are far more complicated yet much less used. In practice, good quantum numbers can be integers or half integers, therefore, we use real numbers to denote them in this module for simplicity. Independent quantum numbers, such as the particle number and the spin z-component, can coexist at the same time. We use type QuantumNumber to represent the complete set of independent ones for a single basis of a Hilbert space, and type QuantumNumbers to represent the whole quantum numbers for the total bases."
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#QuantumNumber-1",
    "page": "Good quantum numbers",
    "title": "QuantumNumber",
    "category": "section",
    "text": "The abstract type for the complete set of independent good quantum numbers for a single basis.Main features include:function fieldnames: get the names of the quantum numbers\nfunction periods: get the periods of the quantum numbers\narithmetic operations: +, -, *, ^, ⊕, ⊗\nhashable: concrete instances can be used as keys for a dict or a set\niterable: concrete instances are iterable over their values\ncomparable: two concrete instances can be comparedIn particular, QuantumNumber <: AbstractNamedVector{Float64}, all features supported by AbstractNamedVector are also available for QuantumNumber. See also AbstractNamedVector.For convenience, 4 kinds of good quantum numbers are predefined in this module, i.e.SQN: for spin z-component reserved systems\nPQN: for particle number reserved systems\nSPQN: for both particle number and spin-z component reserved systems\nZ2QN: for systems with a Z_2 conservation quantum numberUsers who want to define their own Z_N-like quantum numbers must handle the periodicities in the construction function, otherwise, wrong results will be get when arithmetic operations, such as + or -, are involved. It is highly recommended to use the macro @quantumnumber to define your own concrete QuantumNumbers."
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#QuantumNumbers-1",
    "page": "Good quantum numbers",
    "title": "QuantumNumbers",
    "category": "section",
    "text": "The whole quantum numbers for the total bases, which has three forms\'G\' form: the general form, which has no restriction for the contents of the QuantumNumbers\n\'U\' form: the unitary form, which requires no duplicates in the contents of the QuantumNumbers\n\'C\' form: the canonical form, which requires not only no duplicates but also accending-order storage in the contents of the QuantumNumberUsually, \'G\'-formed and \'U\'-formed QuantumNumberses can be transformed to the corresponding \'C\'-formed ones by the sort function.To achieve high efficiency:The contents of a QuantumNumbers are an homogenous array of a certain kind of concrete QuantumNumbers.\nThe quantum numbers are stored in a compressed form similiar to that of a CSC/CSR sparse matrix.Main features include:function eltype: get the concrete type of the quantum numbers it contains\narithmetic operations: +, -, *, ^, ⊗, ⊕\niterable: various iteration supports, including functions such as iterate, keys, values and pairs\n...For a complete summation of its features, please refer to the manual."
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.qnsbruteforce",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.qnsbruteforce",
    "category": "constant",
    "text": "Choice associated with QuantumNumbers, meaning \'by brute force\'.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.qnscontents",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.qnscontents",
    "category": "constant",
    "text": "Choice associated with QuantumNumbers, meaning \'for contents\'.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.qnscounts",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.qnscounts",
    "category": "constant",
    "text": "Choice associated with QuantumNumbers, meaning \'by counts\'.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.qnsexpansion",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.qnsexpansion",
    "category": "constant",
    "text": "Choice associated with QuantumNumbers, meaning \'for expansion\'.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.qnsindices",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.qnsindices",
    "category": "constant",
    "text": "Choice associated with QuantumNumbers, meaning \'for indices\'.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.qnsindptr",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.qnsindptr",
    "category": "constant",
    "text": "Choice associated with QuantumNumbers, meaning \'by indptr\'.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.qnsmontecarlo",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.qnsmontecarlo",
    "category": "constant",
    "text": "Choice associated with quantumnumbers, meaning \'by Monte Carlo\'.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.PQN",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.PQN",
    "category": "type",
    "text": "PQN(N::Real)\n\nThe concrete QuantumNumber of a quantum system with particle number N conserved.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumber",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumber",
    "category": "type",
    "text": "Abstract type for all concrete quantum numbers for a single basis.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumbers",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumbers",
    "category": "type",
    "text": "QuantumNumbers(form::Char,contents::Vector{<:QuantumNumber},counts::Vector{Int},::typeof(qnscounts))\nQuantumNumbers(form::Char,contents::Vector{<:QuantumNumber},indptr::Vector{Int},::typeof(qnsindptr))\n\nThe whole quantum numbers of the total bases of a Hilbert space. The default constructors construct a QuantumNumbers from a vector of concrete quantum numbers and an vector containing their counts or indptr.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumbers",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumbers",
    "category": "type",
    "text": "QuantumNumbers(qn::QuantumNumber,count::Int=1)\n\nConstruct a QuantumNumbers with one unique quantum number which occurs count times.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumbers-Tuple{OrderedCollections.OrderedDict{#s65,Int64} where #s65<:Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumber}",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumbers",
    "category": "method",
    "text": "QuantumNumbers(od::OrderedDict{<:QuantumNumber,Int})\n\nConstruct a QuantumNumbers from an ordered dict containing concrete quantum numbers and their counts.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumbers-Tuple{OrderedCollections.OrderedDict{#s65,UnitRange{Int64}} where #s65<:Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumber}",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumbers",
    "category": "method",
    "text": "QuantumNumbers(od::OrderedDict{<:QuantumNumber,UnitRange{Int}})\n\nConstruct a QuantumNumbers from an ordered dict containing concrete quantum numbers and their slices.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.SPQN",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.SPQN",
    "category": "type",
    "text": "SPQN(N::Real,Sz::Real)\n\nThe concrete QuantumNumber of a quantum system with both particle number N and spin z-component Sz conserved.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.SQN",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.SQN",
    "category": "type",
    "text": "SQN(Sz::Real)\n\nThe concrete QuantumNumber of a quantum system with spin z-component Sz conserved.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.Z2QN",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.Z2QN",
    "category": "type",
    "text": "Z2QN(N::Real)\n\nThe concrete QuantumNumber of a quantum system with a Z₂-like conserved quantity.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.@quantumnumber-Tuple{Any,Any,Any}",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.@quantumnumber",
    "category": "macro",
    "text": "@quantumnumber typename fieldnames fieldperiods\n\nConstruct a concrete QuantumNumber with the type name being typename, fieldnames specified by fieldnames and periods specified by fieldperiods.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.:⊕-Union{Tuple{Tuple{Vararg{#s65,N}} where #s65<:QuantumNumber}, Tuple{N}, Tuple{Tuple{Vararg{#s61,N}} where #s61<:QuantumNumber,Tuple{Vararg{Int64,N}}}} where N",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.:⊕",
    "category": "method",
    "text": "⊕(qns::NTuple{N,<:QuantumNumber},signs::NTuple{N,Int}=ntuple(i->1,N)) where N -> QuantumNumbers\n⊕(qnses::NTuple{N,QuantumNumbers{QN}},signs::NTuple{N,Int}=ntuple(i->1,N)) where {N,QN<:QuantumNumber} -> QuantumNumbers{QN}\n\nGet the direct sum of some QuantumNumbers or QuantumNumberses.\n\nnote: Note\nPhysically, the direct sum of a couple of QuantumNumbers or QuantumNumberses is defined by the direct sum of the bases of the Hilbert spaces they represent. Therefore, the input QuantumNumbers or QuantumNumberses must be homogenous. Inhomogenous \'QuantumNumber\'s must be direct producted first to ensure homogenity before the direct sum.\nApparently, the dimension of the result equals the summation of those of the inputs, which means, even for QuantumNumbers, the result will be naturally a QuantumNumbers because the dimension of the result is largeer than 1.\nSigns of QuantumNumbers or QuantumNumberses can be provided when getting their direct sums.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.:⊗-Union{Tuple{QN}, Tuple{Type{QN},QuantumNumber,QuantumNumber}} where QN<:Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumber",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.:⊗",
    "category": "method",
    "text": "⊗(::Type{QN},qn1::QuantumNumber,qn2::QuantumNumber) where QN<:QuantumNumber -> QN\n⊗(qns::NTuple{N,<:QuantumNumber},signs::NTuple{N,Int}=ntuple(i->1,N)) where N -> QuantumNumber\n⊗(qnses::NTuple{N,QuantumNumbers{QN}},signs::NTuple{N,Int}=ntuple(i->1,N)) where {N,QN<:QuantumNumber} -> QuantumNumbers{QN}\n\nGet the direct product of some QuantumNumbers or QuantumNumberses.\n\nnote: Note\nPhysically, the direct product of a couple of QuantumNumbers or QuantumNumberses are defined by the direct product of the bases of the Hilbert spaces they represent. Therefore, QuantumNumbers with differenct types or QuantumNumberses with differenct eltypes are allowed to be direct producted in principle. However, for simplicity, we only implement a method which handle the situation of two QuantumNumbers with differenct types. The type of the result should be provided as the first parameter. Note that in this situation, the fieldnames and periods of the result type must be exactly equal to the flattened fieldnames and periods of the two input QuantumNumbers, which means, even the order of the input QuantumNumbers matters.\nApparently, the dimension of the result equals the product of those of the inputs. Therefore, the direct product of QuantumNumbers is also a QuantumNumber since its dimension is still one.\nFor other situations except the one mentioned in Note.1, the input QuantumNumbers or QuantumNumberses must be homogenous. Meanwhile, signs can also be provided for these situations. Note that each quantum number in the contents of the result is obtained by a summation of the corresponding quanum numbers out of the inputs with the correct signs. This is a direct observation of the Abelian nature of our quantum numbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.PQNS-Tuple{Real}",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.PQNS",
    "category": "method",
    "text": "PQNS(N::Real)\n\nConstruct the QuantumNumbers of the Hilbert space of a single-particle state with at most N identical particles.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.SPQNS-Tuple{Real}",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.SPQNS",
    "category": "method",
    "text": "SPQNS(S::Real)\n\nConstruct the QuantumNumbers of the Hilbert space of a single site with internal degrees of freedom that can be ascribed to a spin S.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.SQNS-Tuple{Real}",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.SQNS",
    "category": "method",
    "text": "SQNS(S::Real)\n\nConstruct the QuantumNumbers of the Hilbert space of a signle spin S.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.SzPQNS-Tuple{Real}",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.SzPQNS",
    "category": "method",
    "text": "SzPQNS(Sz::Real)\n\nConstruct the QuantumNumbers of the Hilbert space of a single-paritcle state with at most one particle whose spin-z component is Sz.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.Z2QNS-Tuple{}",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.Z2QNS",
    "category": "method",
    "text": "Z2QNS()\n\nConstruct the QuantumNumbers of a Z_2 Hilbert space.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.decompose-Union{Tuple{QN}, Tuple{N}, Tuple{Tuple{Vararg{QuantumNumbers{QN},N}},QN,Tuple{Vararg{Int64,N}},Val{6}}} where QN<:Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumber where N",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.decompose",
    "category": "method",
    "text": "decompose(qnses::NTuple{N,QuantumNumbers{QN}},target::QN,signs::NTuple{N,Int},::typeof(qnsbruteforce);nmax::Int=20) where {N,QN<:QuantumNumber} -> Vector{NTuple{N,Int}}\ndecompose(qnses::NTuple{N,QuantumNumbers{QN}},target::QN,signs::NTuple{N,Int},::typeof(qnsmontecarlo);nmax::Int=20) where {N,QN<:QuantumNumber} -> Vector{NTuple{N,Int}}\n\nFind a couple of decompositions of target with respect to qnses.\n\nnote: Note\nA tuple of integers (i₁,i₂,...) is called a decomposition of a given target with respect to the given qnses if and only if they satisfy the \"decomposition rule\":sum_textj textsignstextjtimestextqnsestextjtexti_textj==texttargetThis equation is in fact a kind of a set of restricted linear Diophantine equations. Indeed, our quantum numbers are always discrete Abelian ones and all instances of a concrete QuantumNumber forms a module over the ring of integers. Therefore, each quantum number can be represented as a integral multiple of the unit element of the Abelian module, which results in the final reduction of the above equation to a set of linear Diophantine equations. Then finding a decomposition is equivalent to find a solution of the reduced linear Diophantine equations, with the restriction that the quantum numbers constructed from the solution should be in the corresponding qnses. Here we provide two methods to find such decompositions, one is by brute force (qnsbruteforce case), and the other is by Monte Carlo simultatioins (qnsmontecarlo case).\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.dimension-Tuple{Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumbers}",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.dimension",
    "category": "method",
    "text": "dimension(qns::QuantumNumbers) -> Int\n\nThe dimension of the Hilbert space a QuantumNumbers represents.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.dimension-Tuple{Type{#s68} where #s68<:Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumber}",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.dimension",
    "category": "method",
    "text": "dimension(::Type{<:QuantumNumber}) -> Int\ndimension(::QuantumNumber) -> Int\n\nThe dimension of the Hilbert space a QuantumNumber represents. Apparently, this is always 1.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.expand-Tuple{Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumbers,Val{3}}",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.expand",
    "category": "method",
    "text": "expand(qns::QuantumNumbers,::typeof(qnscontents)) -> Vector{qns|>eltype}\nexpand(qns::QuantumNumbers,::typeof(qnsindices)) -> Vector{Int}\n\nExpand the contents (qnscontents case) or indices (qnsindices case) of a QuantumNumbers to the uncompressed form.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.regularize!-Union{Tuple{QN}, Tuple{Type{QN},AbstractArray{Float64,1}}} where QN<:Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumber",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.regularize!",
    "category": "method",
    "text": "regularize!(::Type{QN},array::AbstractVector{Float64}) where QN<:QuantumNumber\nregularize!(::Type{QN},array::AbstractMatrix{Float64}) where QN<:QuantumNumber\n\nRegularize the elements of an array in place so that it can represent quantum numbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.regularize-Union{Tuple{QN}, Tuple{Type{QN},Union{AbstractArray{Float64,1}, AbstractArray{Float64,2}}}} where QN<:Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumber",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.regularize",
    "category": "method",
    "text": "regularize(::Type{QN},array::Union{AbstractVector{Float64},AbstractMatrix{Float64}}) where {QN<:QuantumNumber}\n\nRegularize the elements of an array and return a copy that can represent quantum numbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.reorder-Tuple{Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumbers,Array{Int64,1},Val{3}}",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.reorder",
    "category": "method",
    "text": "reorder(qns::QuantumNumbers,permutation::Vector{Int},::typeof(qnscontents)) -> QuantumNumbers\nreorder(qns::QuantumNumbers,permutation::Vector{Int},::typeof(qnsexpansion)) -> QuantumNumbers\n\nReorder the quantum numbers contained in a QuantumNumbers with a permutation and return the new one. For qnscontents case, the permutation is for the contents of the original QuantumNumbers while for qnsexpansion case, the permutation is for the expansion of the original QuantumNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.subset-Union{Tuple{QN}, Tuple{QuantumNumbers{QN},QN}} where QN<:Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumber",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.subset",
    "category": "method",
    "text": "subset(qns::QuantumNumbers{QN},target::QN) where QN<:QuantumNumber -> QuantumNumbers{QN}\nsubset(qns::QuantumNumbers{QN},targets::NTuple{N,QN}) where {N,QN<:QuantumNumber} -> QuantumNumbers{QN}\n\nFind a subset of a QuantumNumbers by picking out the quantum numbers in targets.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.toordereddict-Tuple{Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumbers,Val{1}}",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.toordereddict",
    "category": "method",
    "text": "toordereddict(qns::QuantumNumbers,::typeof(qnsindptr)) -> OrderedDict{qns|>eltype,UnitRange{Int}}\ntoordereddict(qns::QuantumNumbers,::typeof(qnscounts)) -> OrderedDict{qns|>eltype,Int}\n\nConvert a QuantumNumbers to an ordered dict.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.ukron-Union{Tuple{Tuple{Vararg{QuantumNumbers{QN},N}}}, Tuple{QN}, Tuple{N}, Tuple{Tuple{Vararg{QuantumNumbers{QN},N}},Tuple{Vararg{Int64,N}}}} where QN<:Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumber where N",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.ukron",
    "category": "method",
    "text": "ukron(qnses::NTuple{N,QuantumNumbers{QN}},signs::NTuple{N,Int}=ntuple(i->1,N)) where {N,QN<:QuantumNumber} -> QuantumNumbers{QN},Dict{QN,Dict{NTuple{N,QN},UnitRange{Int}}}\n\nUnitary Kronecker product of several QuantumNumberses. The product result as well as the records of the product will be returned.\n\nnote: Note\nAll input QuantumNumbers must be \'U\' formed or \'C\' formed.\nSince duplicate quantum number are not allowed in \'U\' formed and \'C\' formed QuantumNumberses, in general, there exists a merge process of duplicate quantum numbers in the product result. Therefore, records are needed to keep track of this process, which will be returned along with the product result. The records are stored in a Dict{QN,Dict{NTuple{N,QN},UnitRange{Int}}} typed dict, in which, for each unduplicate quantum number qn in the product result, there exist a record Dict((qn₁,qn₂,...)=>start:stop,...) telling what quantum numbers (qn₁,qn₂,...) a mereged duplicate qn comes from and what slice start:stop this merged duplicate corresponds.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.:*-Tuple{Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumber,Integer}",
    "page": "Good quantum numbers",
    "title": "Base.:*",
    "category": "method",
    "text": "*(qn::QuantumNumber,factor::Integer) -> QuantumNumber\n*(factor::Integer,qn::QuantumNumber) -> QuantumNumber\n*(qns::QuantumNumbers,factor::Integer) -> QuantumNumbers\n*(factor::Integer,qns::QuantumNumbers) -> QuantumNumbers\n\nOverloaded * operator for the multiplication between an integer and a QuantumNumber or a QuantumNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.:+-Tuple{Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumber}",
    "page": "Good quantum numbers",
    "title": "Base.:+",
    "category": "method",
    "text": "+(qn::QuantumNumber) -> QuantumNumber\n+(qn::QN,qns::QN...) where QN<:QuantumNumber -> QN\n+(qns::QuantumNumbers) -> QuantumNumbers\n+(qn::QN,qns::QuantumNumbers{QN}) where QN<:QuantumNumber -> QuantumNumbers{QN}\n+(qns::QuantumNumbers{QN},qn::QN) where QN<:QuantumNumber -> QuantumNumbers{QN}\n\nOverloaded + operator for QuantumNumber and QuantumNumbers.\n\nnote: Note\nThe addition between a QuantumNumbers and a QuantumNumber is just a global shift of the contents of the QuantumNumbers by the QuantumNumber, therefore, the result is a QuantumNumbers.\n+ cannot be used between two QuantumNumbers because the result is ambiguous. Instead, use ⊕ for direct sum and ⊗ for direct product.\nTo ensure type stability, two QuantumNumber can be added together if and only if they are of the same type.\nSimilarly, a QuantumNumber and a QuantumNumbers can be added together if and only if the former\'s type is the same with the latter\'s eltype.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.:--Tuple{Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumber}",
    "page": "Good quantum numbers",
    "title": "Base.:-",
    "category": "method",
    "text": "-(qn::QuantumNumber) -> QuantumNumber\n-(qn1::QN,qn2::QN) where QN<:QuantumNumber -> QN\n-(qns::QuantumNumbers) -> QuantumNumbers\n-(qn::QN,qns::QuantumNumbers{QN}) where QN<:QuantumNumber -> QuantumNumbers{QN}\n-(qns::QuantumNumbers{QN},qn::QN) where QN<:QuantumNumber -> QuantumNumbers{QN}\n\nOverloaded - operator for QuantumNumber and QuantumNumbers.\n\nnote: Note\nThe subtraction between a QuantumNumbers and a QuantumNumber is just a global shift of the contents of the QuantumNumbers by the QuantumNumber, therefore, the result is a QuantumNumbers.\n- cannot be used between two QuantumNumbers because the result is ambiguous. Instead, use ⊕ with signs for direct sum and ⊗ with signs for direct product.\nTo ensure type stability, a QuantumNumber can be subtracted by another QuantumNumber if and only if they are of the same type.\nSimilarly, a QuantumNumber can be subtracted by a QuantumNumbers or vice versa if and only if the former\'s type is the same with the latter\'s eltype.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.:==-Tuple{Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumbers,Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumbers}",
    "page": "Good quantum numbers",
    "title": "Base.:==",
    "category": "method",
    "text": "==(qns1::QuantumNumbers,qns2::QuantumNumbers) -> Bool\n\nOverloaded == operator. Two QuantumNumberses are equal to each other if and only if both their contentses and indptrs are elementwise equal to each other.\n\nnote: Note\nIt is not necessary for two QuantumNumberses to have the same eltype nor the same form to be equal to each other.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.:^-Tuple{Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumber,Integer}",
    "page": "Good quantum numbers",
    "title": "Base.:^",
    "category": "method",
    "text": "^(qn::QuantumNumber,factor::Integer) -> QuantumNumber\n^(qns::QuantumNumbers,factor::Integer) -> QuantumNumbers\n\nOverloaded ^ operator for QuantumNumber and QuantumNumbers. This operation translates into the direct product ⊗ of factor copies of qn or qns.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.eltype-Union{Tuple{Type{#s68} where #s68<:QuantumNumbers{QN}}, Tuple{QN}} where QN",
    "page": "Good quantum numbers",
    "title": "Base.eltype",
    "category": "method",
    "text": "eltype(::Type{<:QuantumNumbers{QN}}) where QN\neltype(qns::QuantumNumbers)\n\nGet the type of the concrete QuantumNumber contained in a QuantumNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.findall-Union{Tuple{QN}, Tuple{QuantumNumbers{QN},QN,Union{Val{3}, Val{4}}}} where QN<:Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumber",
    "page": "Good quantum numbers",
    "title": "Base.findall",
    "category": "method",
    "text": "findall(qns::QuantumNumbers{QN},target::QN,choice::Union{typeof(qnscontents),typeof(qnsexpansion)}) where QN<:QuantumNumber -> Vector{Int}\nfindall(qns::QuantumNumbers{QN},targets::NTuple{N,QN},::typeof(qnscontents)) where {N,QN<:QuantumNumber} -> Vector{Int}\nfindall(qns::QuantumNumbers{QN},targets::NTuple{N,QN},::typeof(qnsexpansion)) where {N,QN<:QuantumNumber} -> Vector{Int}\n\nFind all the indices of the target quantum numbers in the contents (qnscontents case) or the expansion (qnsexpansion case) of a QuantumNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.getindex-Tuple{Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumbers,Int64}",
    "page": "Good quantum numbers",
    "title": "Base.getindex",
    "category": "method",
    "text": "getindex(qns::QuantumNumbers,index::Int) -> QuantumNumber\ngetindex(qns::QuantumNumbers,slice::UnitRange{Int}) -> QuantumNumbers\ngetindex(qns::QuantumNumbers,indices::Vector{Int}) -> QuantumNumbers\n\nOverloaded [] operator.\n\nnote: Note\nFor a QuantumNumbers, all these getindex functions act on its contents, i.e. its compressed data, but not on its expansion, i.e. the uncompressed data. This definition is consistent with the length function.\nWhen the index is an integer, the result is a QuantumNumber, while when the index is a unit range or a vector of intgers, the result is a QuantumNumbers. The logic is quite reasonable because such behaviors are much alike to those of a vector container.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.iterate",
    "page": "Good quantum numbers",
    "title": "Base.iterate",
    "category": "function",
    "text": "iterate(qns::QuantumNumbers,state::Int=1)\niterate(rv::Iterators.Reverse{<:QuantumNumbers},state::Int=length(rv.itr,false))\n\nIterate or reversely iterate over the concrete QuantumNumbers contained in a QuantumNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.keys-Tuple{Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumbers}",
    "page": "Good quantum numbers",
    "title": "Base.keys",
    "category": "method",
    "text": "keys(qns::QuantumNumbers) -> Vector{qns|>eltype}\n\nIterate over the concrete QuantumNumbers contained in a QuantumNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.kron-Union{Tuple{QN}, Tuple{Type{QN},QuantumNumber,QuantumNumber}} where QN<:Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumber",
    "page": "Good quantum numbers",
    "title": "Base.kron",
    "category": "method",
    "text": "kron(::Type{QN},qn1::QuantumNumber,qn2::QuantumNumber) where QN<:QuantumNumber -> QN\nkron(qns::NTuple{N,<:QuantumNumber},signs::NTuple{N,Int}=ntuple(i->1,N)) where N -> QuantumNumber\nkron(qnses::NTuple{N,QuantumNumbers{QN}},signs::NTuple{N,Int}=ntuple(i->1,N)) where {N,QN<:QuantumNumber} -> QuantumNumbers{QN}\n\nKronecker product of some QuantumNumbers or QuantumNumberses. This is defined to be equivalent to the direct product ⊗.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.length-Tuple{Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumbers}",
    "page": "Good quantum numbers",
    "title": "Base.length",
    "category": "method",
    "text": "length(qns::QuantumNumbers) -> Int\n\nGet the number of unduplicate qunatum numbers in the QuantumNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.pairs-Tuple{Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumbers,Union{Val{1}, Val{2}}}",
    "page": "Good quantum numbers",
    "title": "Base.pairs",
    "category": "method",
    "text": "pairs(qns::QuantumNumbers,choice::Union{typeof(qnsindptr),typeof(qnscounts)})\n\nIterate over the QuantumNumber=>slice or QuantumNumber=>count pairs.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.show-Tuple{IO,Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumbers}",
    "page": "Good quantum numbers",
    "title": "Base.show",
    "category": "method",
    "text": "show(io::IO,qns::QuantumNumbers)\n\nShow a QuantumNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.sort-Tuple{Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumbers}",
    "page": "Good quantum numbers",
    "title": "Base.sort",
    "category": "method",
    "text": "sort(qns::QuantumNumbers) -> QuantumNumbers,Vector{Int}\n\nSort the quantum numbers of a QuantumNumber, return the sorted QuantumNumber and the permutation array that sorts the expansion of the original QuantumNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.string-Tuple{Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumbers}",
    "page": "Good quantum numbers",
    "title": "Base.string",
    "category": "method",
    "text": "string(qns::QuantumNumbers) -> String\n\nConvert a QuantumNumbers to string.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.values-Tuple{Hamiltonian.Utilities.GoodQuantumNumber.QuantumNumbers,Val{1}}",
    "page": "Good quantum numbers",
    "title": "Base.values",
    "category": "method",
    "text": "values(qns::QuantumNumbers,::typeof(qnsindptr))\nvalues(qns::QuantumNumbers,::typeof(qnscounts))\n\nIterate over the slices/counts of the QuantumNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#qnmanual-1",
    "page": "Good quantum numbers",
    "title": "Manual",
    "category": "section",
    "text": "Modules=[GoodQuantumNumber]\nOrder=  [:module,:constant,:type,:macro,:function]"
},

]}
