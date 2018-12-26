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
    "text": "This module contains the utilities of the Hamiltonian package, all of whose variables will NOT be exported by Hamiltonian."
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
    "location": "man/Utilities/Introduction.html#Hamiltonian.Utilities.forder",
    "page": "Introduction",
    "title": "Hamiltonian.Utilities.forder",
    "category": "constant",
    "text": "forder\n\nIndicate that the convertion between Cartesian index and linear index is using the Fortran order.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Introduction.html#Hamiltonian.Utilities.corder",
    "page": "Introduction",
    "title": "Hamiltonian.Utilities.corder",
    "category": "constant",
    "text": "corder\n\nIndicate that the convertion between Cartesian index and linear index is using the C/C++ order.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Introduction.html#Hamiltonian.Utilities.ind2sub",
    "page": "Introduction",
    "title": "Hamiltonian.Utilities.ind2sub",
    "category": "function",
    "text": "ind2sub(dims::Tuple,ind::Int,order::FOrder) -> Tuple\nind2sub(dims::Tuple,ind::Int,order::COrder) -> Tuple\n\nConvert an linear index to Cartesian index. Fortran-order or C-order can be assigned.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Introduction.html#Hamiltonian.Utilities.sub2ind",
    "page": "Introduction",
    "title": "Hamiltonian.Utilities.sub2ind",
    "category": "function",
    "text": "sub2ind(dims::NTuple{N,Int},inds::NTuple{N,Int},order::FOrder) where N -> Int\nsub2ind(dims::NTuple{N,Int},inds::NTuple{N,Int},order::COrder) where N -> Int\n\nConvert an Cartesian index to linear index. Fortran-order or C-order can be assigned.\n\n\n\n\n\n"
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
    "location": "man/Utilities/Introduction.html#Hamiltonian.Utilities.efficientoperations",
    "page": "Introduction",
    "title": "Hamiltonian.Utilities.efficientoperations",
    "category": "constant",
    "text": "efficientoperations\n\nIndicate that the efficient operations, i.e. \"==\"/\"isequal\", \"<\"/\"isless\" or \"replace\", will be used.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Introduction.html#Hamiltonian.Utilities.delta",
    "page": "Introduction",
    "title": "Hamiltonian.Utilities.delta",
    "category": "function",
    "text": "delta(i,j) -> Int\n\nKronecker delta function.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Introduction.html#Useful-constants-and-functions-1",
    "page": "Introduction",
    "title": "Useful constants and functions",
    "category": "section",
    "text": "atol\nrtol\nFloat\nforder\ncorder\nind2sub\nsub2ind\ndecimaltostr\nordinal\nefficientoperations\ndelta"
},

{
    "location": "man/Utilities/Introduction.html#Hamiltonian.Utilities.:⊕",
    "page": "Introduction",
    "title": "Hamiltonian.Utilities.:⊕",
    "category": "function",
    "text": "Generic interface of the direct sum of some types.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Introduction.html#Hamiltonian.Utilities.:⊗",
    "page": "Introduction",
    "title": "Hamiltonian.Utilities.:⊗",
    "category": "function",
    "text": "Generic interface of the direct product of some types.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Introduction.html#Hamiltonian.Utilities.rank",
    "page": "Introduction",
    "title": "Hamiltonian.Utilities.rank",
    "category": "function",
    "text": "Generic interface of the rank of some types.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Introduction.html#Hamiltonian.Utilities.dimension",
    "page": "Introduction",
    "title": "Hamiltonian.Utilities.dimension",
    "category": "function",
    "text": "Generic interface of the dimension of some types.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Introduction.html#Hamiltonian.Utilities.expand",
    "page": "Introduction",
    "title": "Hamiltonian.Utilities.expand",
    "category": "function",
    "text": "Generic interface of the expansion of some types.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Introduction.html#Hamiltonian.Utilities.permute",
    "page": "Introduction",
    "title": "Hamiltonian.Utilities.permute",
    "category": "function",
    "text": "Generic interface of the permutation of some types.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Introduction.html#Hamiltonian.Utilities.vector",
    "page": "Introduction",
    "title": "Hamiltonian.Utilities.vector",
    "category": "function",
    "text": "Generic interface of the vector representation of some types.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Introduction.html#Hamiltonian.Utilities.matrix",
    "page": "Introduction",
    "title": "Hamiltonian.Utilities.matrix",
    "category": "function",
    "text": "Generic interface of the matrix representation of some types.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Introduction.html#Generic-functions-for-overloading-1",
    "page": "Introduction",
    "title": "Generic functions for overloading",
    "category": "section",
    "text": "⊕\n⊗\nrank\ndimension\nexpand\npermute\nvector\nmatrix"
},

{
    "location": "man/Utilities/Introduction.html#Prerequisites-for-Essentials-1",
    "page": "Introduction",
    "title": "Prerequisites for Essentials",
    "category": "section",
    "text": "Pages=[\n    \"Factory.md\",\n    \"CompositeStructure.md\",\n    \"Tree.md\",\n    \"NamedVector.md\",\n    \"AlgebraOverField.md\",\n    ]\nDepth=2"
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
    "text": "The aim of Factory is to provide tools to hack into Julia codes without knowing the details of their abstract syntax trees and regularize the mechanism to \"escape\" variables in Expr expressions, so that users can manipulate the existing codes, modify them and generate new ones in macros. In particular, Factory in this module means the representation of certain blocks of Julia codes by a usual Julia struct. This representation is much easier to comprehend than the canonical Expr representation and makes it far more convenient to define macros. In general, we propose the following requirements that any factory must satisfy:DECOMPOSITION - An Expr expression can be decomposed into its corresponding factory by the factory\'s constructor.\nCOMPOSITION - A factory can compose its corresponding Expr expression by calling itself.\nESCAPE - A variable should be or not be escaped in the composed Expr expression by a factory depends on predefined escape mechanisms.These three requirements also define the basic interfaces to interact with factories. In practice, we combine the second and third in a single interface, i.e. by passing an instance of certain concrete EscapeMechanism as the only argument of calling a factory, the needed Expr expression with variables correctly escaped can be obtained."
},

{
    "location": "man/Utilities/Factory.html#Escape-mechanisms-1",
    "page": "Factory",
    "title": "Escape mechanisms",
    "category": "section",
    "text": "We adopt Julia structs to denote escape mechanisms so that we can utilize Julia\'s multidispatch to implement different mechanisms whereas keeping the same interface."
},

{
    "location": "man/Utilities/Factory.html#EscapeMechanism-1",
    "page": "Factory",
    "title": "EscapeMechanism",
    "category": "section",
    "text": "EscapeMechanism is the abstract type for all concrete escape mechanisms."
},

{
    "location": "man/Utilities/Factory.html#Escaped-1",
    "page": "Factory",
    "title": "Escaped",
    "category": "section",
    "text": "Escaped has only one attribute:names::NTuple{N,Symbol} where N: the names of variables to be escapedApprently, a variable should be escaped if its name is in the names of an Escaped. This mechanism suits a factory whose variables should be unescaped by default."
},

{
    "location": "man/Utilities/Factory.html#UnEscaped-1",
    "page": "Factory",
    "title": "UnEscaped",
    "category": "section",
    "text": "UnEscaped also has only on attribute:names::NTuple{N,Symbol} where N: the names of variables not to be escapedObviously, on the contrary to [Escaped], a variable should be escaped if its name is not in the names of an UnEscaped. This mechanism suits a factory whose variables should be escaped by default."
},

{
    "location": "man/Utilities/Factory.html#MixEscaped-1",
    "page": "Factory",
    "title": "MixEscaped",
    "category": "section",
    "text": "MixEscaped has two attributes:escaped::Escaped: the escaped part of the mixed mechanism\nunescaped::UnEscaped: the UnEscaped part of the mixed mechanismThis mechanism suits complex factories that parts of it suit the \"escaped\" mechanism while others suit the \"unescaped\" mechanism."
},

{
    "location": "man/Utilities/Factory.html#RawExpr-1",
    "page": "Factory",
    "title": "RawExpr",
    "category": "section",
    "text": "RawExpr has no attributes and it means \"raw expression without any variable escaped\". This mechanism is used for the print of all factories by default."
},

{
    "location": "man/Utilities/Factory.html#Concrete-factories-1",
    "page": "Factory",
    "title": "Concrete factories",
    "category": "section",
    "text": "Out of practical purposes, we implemente 7 kinds of factories, i.e. Inference, Argument, Parameter, Field, Block, FunctionFactory and TypeFactory, which represent a type inference, a function argument, a method or type parameter, a struct field, a begin ... end block, a function itself and a struct itself, respectively. Some of the basic methods making the above requirements fulfilled with these types are based on the powerful functions defined in MacroTools.We want to give a remark that although the types and functions provided in this module helps a lot for the definition of macros, macros should not be abused. On the one hand, some macros may change the language specifications, which makes it hard to understand the codes, and even splits the community; on the one hand, macros usually increases the precompiling/jit time, which means enormous uses of macros in a module may lead to an extremely long load time. Besides, due to the limited ability of the author, the codes in this module are not optimal, which adds to the jit overhead. Any promotion that keeps the interfaces unchanged is welcomed."
},

{
    "location": "man/Utilities/Factory.html#Inference-1",
    "page": "Factory",
    "title": "Inference",
    "category": "section",
    "text": "An Inference has 3 attributes:head::Union{Symbol,Nothing}: the head of the type inference, which must be one of (nothing,:(::),:(<:),:curly)\nname::Union{Symbol,Nothing}: the name of the type inference\nparams::Union{Inference,Vector{Inference},Nothing}: the parameters of the type inferenceAll valid expressions representing type inferences can be passed to the constructor:Inference(:T)\nInference(:(<:Number))\nInference(:(Vector{T}))\nInference(:(Vector{Tuple{String,Int}}))\nInference(:(Type{<:Number}))On the other hand, you can use the macro @inference to construct an Inference directly from a type inference:@inference Vector{Tuple{String,Int}}note: Note\nInference is a recursive struct, i.e. it recursively decomposes a type inference until the final type inference is just a Symbol.\nWhen the input expression is a Symbol, the head and params attributes of the resulting Inference is nothing. Otherwise, its head is the same with that of the input expression, and the args of the input expression will be further decomposed, whose result will be stored in params.\nWhen the head of the input expression is :(<:), the params is an Inference whereas when the head of the input expression is :curly, the params is a Vector{Inference}.Inference uses the UnEscaped mechanism to escape variables, e.g.Inference(:(Vector{T}))(UnEscaped()) |> println\nInference(:(Vector{T}))(UnEscaped(:T)) |> println\nInference(:(Vector{T}))(UnEscaped(:Vector,:T)) |> println"
},

{
    "location": "man/Utilities/Factory.html#Argument-1",
    "page": "Factory",
    "title": "Argument",
    "category": "section",
    "text": "An Argument has 4 attributes:name::Union{Symbol,Nothing}: the name of the argument\ntype::Inference: the type inference of the argument\nslurp::Bool: whether the argument should be expanded by ...\ndefault::Any: the default value of the argument, nothing for those with no default valuesAll valid expressions representing the arguments of functions can be passed to the constructor:Argument(:arg)\nArgument(:(arg::ArgType))\nArgument(:(arg::ArgType...))\nArgument(:(arg::ArgType=default))Or you can use the macro @argument for a direct construction from an argument declaration:@argument arg::ArgType=defaultThe construction from such expressions is based on the the MacroTools.splitarg function.Argument uses the MixEscaped mechanism to escape variables, with the UnEscaped mechanism for type and Escaped mechanism for default, e.g.Argument(:(arg::Real=zero(Int)))(MixEscaped(UnEscaped(),Escaped(:zero,:Int))) |> printlnIt can be seen the name of an argument will never be escaped, which is obvious since the name of a function argument is always local. By the way, the composition of an Argument expression is based on the MacroTools.combinearg function."
},

{
    "location": "man/Utilities/Factory.html#Parameter-1",
    "page": "Factory",
    "title": "Parameter",
    "category": "section",
    "text": "A Parameter has 2 attributes:name::Union{Symbol,Nothing}: the name of the parameter\ntype::Union{Inference,Nothing}: the type inference of the parameterAll expressions that represent type parameters or method parameters are allowed to be passed to the constructor:Parameter(:T)\nParameter(:(<:Number))\nParameter(:(T<:Number))The macro @parameter completes the construction directly from a parameter declaration:@parameter T<:Numbernote: Note\nWe use nothing to denote a missing name or type.\nTwo subtle situations of type/method parameters, e.g. MyType{T} and MyType{Int}, should be distinguished by Parameter(:T) and Parameter(:(<:Int)). In other words, MyType{Int} is in fact not supported. Indeed, Parameter(:Int) will treat :Int as the parameter name but not the parameter type.Parameter uses the UnEscaped mechanism to escape variables, too, e.g.Parameter(:(N<:Vector{T}))(UnEscaped(:T)) |> printlnAs is similar to Argument, the name of a method/type parameter will never be escaped because of its local scope."
},

{
    "location": "man/Utilities/Factory.html#Field-1",
    "page": "Factory",
    "title": "Field",
    "category": "section",
    "text": "A Field has 2 attributes:name::Symbol: the name of the field\ntype::Inference: the type inference of the fieldLegal expressions can be used to construct a Field instance by its constructor:Field(:field)\nField(:(field::FieldType))\nField(:(field::ParametricType{T}))The macro @field is also provided to help the construction directly from a field declaration:@field field::FieldTypeThe construction from these expressions is based on the MacroTools.splitarg function.Field uses the UnEscaped mechanism to escape variables as well, e.g.Field(:(field::Dict{N,D}))(UnEscaped(:N,:D)) |> printlnThe name of a struct will never be escaped either because it is a local variable tightly binding to a struct. It is noted that the composition of field expressions is based on the MacroTools.combinefield function."
},

{
    "location": "man/Utilities/Factory.html#Block-1",
    "page": "Factory",
    "title": "Block",
    "category": "section",
    "text": "A Block has only one attribute:body::Vector{Any}: the body of the begin ... end blockAny expression can be passed to the constructor of Block:Block(:(x=1))\nBlock(:(x=1;y=2))\nBlock(:(begin x=1 end))\nBlock(quote\n        x=1\n        y=2\n    end)Or you can construct a Block instance directly from any code by the macro @block:@block x=1 y=2The body of a block can also be extended by the push! function or the @push! macro.note: Note\nThe body of a Block is somewhat \"flattened\", i.e. it contains no begin ... end blocks. During the initialization, any such input block will be unblocked and added to the body part by part. So is the push! and @push! procedures.\nAll LineNumberNodes generated by the input codes will also be included in the block\'s body. However, you can use rmlines! or @rmlines! to remove them from the body of an existing Block, or use rmlines or @rmlines to get a copy with them removed in the body.Block uses the Escaped mechanism to escape variables. This is because variables in a block are often local ones and should not be escaped. Therefore, only those defined in other modules should be noted and escaped, which usually constitute the minority. For example,Block(:(x=1;y=2;z=Int[1,2,3]))(Escaped(:Int)) |> println"
},

{
    "location": "man/Utilities/Factory.html#FunctionFactory-1",
    "page": "Factory",
    "title": "FunctionFactory",
    "category": "section",
    "text": "A FunctionFactory has 7 attributes:name::Union{Symbol,Expr}: the name of the function\nparams::Vector{Inference}: the method parameters of the function\nargs::Vector{Argument}: the positional arguments of the function\nkwargs::Vector{Argument}: the keyword arguments of the function\nrtype::Inference: the return type of the function\nwhereparams::Vector{Parameter}: the method parameters specified by the where keyword\nbody::Block: the body of the functionAll expressions that represent functions are allowed to be passed to the constructor:FunctionFactory(:(f()=nothing))\nFunctionFactory(:(f(x)=x))\nFunctionFactory(:(f(x::Int,y::Int;choice::Function=sum)=choice(x,y)))\nFunctionFactory(:(f(x::T,y::T;choice::Function=sum) where T<:Number=choice(x,y)))\nFunctionFactory(:((f(x::T,y::T;choice::Function=sum)::T) where T<:Number=choice(x,y)))\nFunctionFactory(:(\n    function (f(x::T,y::T;choice::Function=sum)::T) where T<:Number\n        choice(x,y)\n    end\n))\nFunctionFactory(\n    quote\n        function (f(x::T,y::T;choice::Function=sum)::T) where T<:Number\n            choice(x,y)\n        end\n    end\n)Similarly, an instance can also be constructed from the macro @functionfactory:@functionfactory (f(x::T,y::T;choice::Function=sum)::T) where T<:Number=choice(x,y)The construction from such expressions are based on the MacroTools.splitdef function.note: Note\nSince Julia 0.7, the form MyType{D}(data::D) where D only appears in struct constructors, therefore, the attribute :params of a function factory is nonempty only when this factory aims to represent a struct constructor.\nUsually, the name of a function factory is a Symbol. However, if the factory aims to extend some methods of a function defined in another module, e.g., Base.eltype, the name will be an Expr.Since FunctionFactory adopts the MixEscaped mechanism to escape variables, with UnEscaped for params, args, kwargs, rtype and whereparams while Escaped for name and body. It is worth to emphasize that the name of a function factory belongs to the Escaped part. Therefore, when it is an Expr, it will never be escaped because an Expr cannot be a element of a NTuple{N,Symbol} where N. See examples,FunctionFactory(:((f(x::T,y::T;choice::Function=sum)::T) where T<:Number=max(x,y,choice(x,y))))(MixEscaped(UnEscaped(:T),Escaped(:f,:max,))) |> println\nFunctionFactory(:((f(x::T,y::T;choice::Function=sum)::T) where T<:Number=max(x,y,choice(x,y))))(MixEscaped(UnEscaped(:T),Escaped(:max))) |> printlnThe compositions of function expressions are based on the MacroTools.combinedef function.Other features include:Positional arguments can be added by addargs! or @addargs!\nKeyword arguments can be added by addkwargs! or @addkwargs!\nWhere parameters can be added by addwhereparams! or @addwhereparams!\nBody can be extended by extendbody! or @extendbody!"
},

{
    "location": "man/Utilities/Factory.html#TypeFactory-1",
    "page": "Factory",
    "title": "TypeFactory",
    "category": "section",
    "text": "A TypeFactory has 6 attributes:name::Symbol: the name of the struct\nmutable::Bool: whether or not the struct is mutable\nparams::Vector{Parameter}: the type parameters of the struct\nsupertype::Inference: the supertype of the struct\nfields::Vector{Field}: the fields of the struct\nconstructors::Vector{FunctionFactory}: the inner constructors of the structAny expression representing valid struct definitions can be passed to the constructor:TypeFactory(:(struct StructName end))\nTypeFactory(:(struct StructName{T} end))\nTypeFactory(:(struct Child{T} <: Parent{T} end))\nTypeFactory(:(\n    struct Child{T<:Number} <: Parent{T}\n        field1::T\n        field2::T\n    end\n))\nTypeFactory(\n    quote\n        struct Child{T<:Number} <: Parent{T}\n            field1::T\n            field2::T\n        end\n    end\n)Also, the macro @typefactory supports the construction directly from a type definition:@typefactory struct Child{T<:Number} <: Parent{T}\n                field1::T\n                field2::T\n                Child(field1::T,field2::T=zero(T)) where T=new{T}(field1,field2)\n            endThe construction from these expressions is based on the MacroTools.splitstructdef function.TypeFactory also uses the MixEscaped mechanism to escape variables, with the UnEscaped part for params, supertype and fields, the Escaped part for name, and both for constructors. For example,@typefactory(struct Child{T<:Number} <: Parent{T}\n    field::T\n    Child(field::T) where T=new{T}(field)\nend)(MixEscaped(UnEscaped(:T),Escaped(:Child))) |>printlnThe composition of a type expression is based on the MacroTools.combinestructdef function.Other features include:Fields can be added by addfields! or @addfields!\nType parameters can be added by addparams! or @addparams!\nInner constructors can be added by addconstructors! or @addconstructors!"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.FExpr",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.FExpr",
    "category": "constant",
    "text": "Factory expression types, which is defined as Union{Symbol,Expr}.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.rawexpr",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.rawexpr",
    "category": "constant",
    "text": "rawexpr\n\nIndicate that no variable in a factory should be escaped.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.Argument",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.Argument",
    "category": "type",
    "text": "Argument(name::Union{Symbol,Nothing},type::Inference,slurp::Bool,default::Any)\nArgument(;name::Union{Symbol,Nothing}=nothing,type::Inference=Inference(:Any),slurp::Bool=false,default::Any=nothing)\nArgument(expr::FExpr)\n\nThe struct to describe a argument of a function.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.Argument-Tuple{Hamiltonian.Utilities.Factory.RawExpr}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.Argument",
    "category": "method",
    "text": "(a::Argument)(em::RawExpr) -> Expr\n(a::Argument)(em::MixEscaped) -> Expr\n\nConvert an Argument to the Expr representation of the argument it describes.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.Block",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.Block",
    "category": "type",
    "text": "Block(parts::FExpr...)\n\nThe struct to describe a begin ... end block.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.Block-Tuple{Hamiltonian.Utilities.Factory.MixEscaped}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.Block",
    "category": "method",
    "text": "(b::Block)(em::RawExpr) -> Expr\n(b::Block)(em::Escaped) -> Expr\n(b::Block)(em::MixEscaped) -> Expr\n\nConvert a Block to the Expr representation of the begin ... end block it describes.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.Escaped",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.Escaped",
    "category": "type",
    "text": "Escaped(names::Symbol...)\n\nIndicate that symbols of a factory should be escaped if they are in names.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.Field",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.Field",
    "category": "type",
    "text": "Field(name::Symbol,type::Inference)\nField(;name::Symbol,type::FExpr=Inference(:Any))\nField(expr::FExpr)\n\nThe struct to describe a field of a struct.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.Field-Tuple{Union{RawExpr, #s14, #s13} where #s13<:Hamiltonian.Utilities.Factory.MixEscaped where #s14<:Hamiltonian.Utilities.Factory.UnEscaped}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.Field",
    "category": "method",
    "text": "(f::Field)(em::RawExpr) -> Expr\n(f::Field)(em::UnEscaped) -> Expr\n(f::Field)(em::MixEscaped) -> Expr\n\nConvert a Field to the Expr representation of the field it describes.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.FunctionFactory",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.FunctionFactory",
    "category": "type",
    "text": "FunctionFactory(name::FExpr,params::Vector{Inference},args::Vector{Argument},kwargs::Vector{Argument},rtype::Inference,whereparams::Vector{Parameter},body::Block)\nFunctionFactory(    ;name::FExpr,\n                    params::Vector{Inference}=Inference[],\n                    args::Vector{Argument}=Argument[],\n                    kwargs::Vector{Argument}=Argument[],\n                    rtype::Inference=Inference(:Any),\n                    whereparams::Vector{Parameter}=Parameter[],\n                    body::Block=Block()\n                    )\nFunctionFactory(expr::Expr)\n\nThe struct to describe a function.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.FunctionFactory-Tuple{Union{RawExpr, #s68} where #s68<:Hamiltonian.Utilities.Factory.MixEscaped}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.FunctionFactory",
    "category": "method",
    "text": "(ff::FunctionFactory)(em::RawExpr) -> Expr\n(ff::FunctionFactory)(em::MixEscaped) -> Expr\n\nConvert a FunctionFactory to the Expr representation of the function it describes.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.Inference",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.Inference",
    "category": "type",
    "text": "Inference(head::Union{Symbol,Nothing},name::Union{Symbol,Nothing},params::Union{Inference,Vector{Inference},Nothing})\nInference(;\n        head::Union{Symbol,Nothing}=nothing,\n        name::Union{Symbol,Nothing}=nothing,\n        params::Union{Inference,Vector{Inference},Nothing}=nothing,\n        )\nInference(expr::FExpr)\n\nThe struct to describe a type inference.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.Inference-Tuple{Hamiltonian.Utilities.Factory.MixEscaped}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.Inference",
    "category": "method",
    "text": "(i::Inference)(em::RawExpr) -> FExpr\n(i::Inference)(em::UnEscaped) -> FExpr\n(i::Inference)(em::MixEscaped) -> FExpr\n\nConvert a Inference to the Expr representation of the type inference it describes.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.MixEscaped",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.MixEscaped",
    "category": "type",
    "text": "MixEscaped(escaped::Escaped)\nMixEscaped(unescaped::UnEscaped)\nMixEscaped(escaped::Escaped,unescaped::UnEscaped)\nMixEscaped(unescaped::UnEscaped,escaped::Escaped)\n\nIndicate that some parts of a factory use the Escaped mechanism while other parts use the UnEscaped mechanism.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.Parameter",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.Parameter",
    "category": "type",
    "text": "Parameter(name::Union{Symbol,Nothing},type::Union{Inference,Nothing})\nParameter(;name::Union{Symbol,Nothing}=nothing,type::Union{Inference,Nothing}=nothing)\nParameter(expr::FExpr)\n\nThe struct to describe a parameter of a function or a type.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.Parameter-Tuple{Union{RawExpr, #s14, #s13} where #s13<:Hamiltonian.Utilities.Factory.MixEscaped where #s14<:Hamiltonian.Utilities.Factory.UnEscaped}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.Parameter",
    "category": "method",
    "text": "(p::Parameter)(em::RawExpr) -> FExpr\n(p::Parameter)(em::UnEscaped) -> FExpr\n(p::Parameter)(em::MixEscaped) -> FExpr\n\nConvert a Parameter to the Expr representation of the parameter it describes.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.RawExpr",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.RawExpr",
    "category": "type",
    "text": "Raw expression without any variable escaped.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.TypeFactory",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.TypeFactory",
    "category": "type",
    "text": "TypeFactory(name::Symbol,mutable::Bool,params::Vector{Parameter},supertype::Inference,fields::Vector{Field},constructors::Vector{FunctionFactory})\nTypeFactory(    ;name::Symbol,\n                mutable::Bool=false,\n                params::Vector{Parameter}=Parameter[],\n                supertype::Inference=Inference(:Any),\n                fields::Vector{Field}=Field[],\n                constructors::Vector{FunctionFactory}=FunctionFactory[],\n                )\nTypeFactory(expr::Expr)\n\nThe struct to describe a struct.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.TypeFactory-Tuple{Union{RawExpr, #s68} where #s68<:Hamiltonian.Utilities.Factory.MixEscaped}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.TypeFactory",
    "category": "method",
    "text": "(tf::TypeFactory)(em::RawExpr) -> Expr\n(tf::TypeFactory)(em::MixEscaped) -> Expr\n\nConvert a TypeFactory to the Expr representation of the struct it describes.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.UnEscaped",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.UnEscaped",
    "category": "type",
    "text": "UnEscaped(names::Symbol...)\n\nIIndicate that symbols of a factory should be escaped if they are not in names.\n\n\n\n\n\n"
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
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.addparams!-Tuple{Hamiltonian.Utilities.Factory.FunctionFactory}",
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
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.escape-Tuple{Any,Hamiltonian.Utilities.Factory.RawExpr}",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.escape",
    "category": "method",
    "text": "escape(expr,::RawExpr) -> Any\nescape(expr,::Escaped) -> Any\nescape(expr,::UnEscaped) -> Any\nescape(expr::Symbol,em::Escaped) -> FExpr\nescape(expr::Expr,em::Escaped) -> Expr\nescape(expr::Symbol,em::UnEscaped) -> FExpr\nescape(expr::Expr,em::UnEscaped) -> Expr\n\nEscape the variables in the input expression.\n\n\n\n\n\n"
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
    "location": "man/Utilities/Factory.html#Hamiltonian.Utilities.Factory.EscapeMechanism",
    "page": "Factory",
    "title": "Hamiltonian.Utilities.Factory.EscapeMechanism",
    "category": "type",
    "text": "Abstract escape mechanism.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Factory.html#Base.:==-Union{Tuple{F}, Tuple{F,F}} where F<:Hamiltonian.Utilities.Factory.AbstractFactory",
    "page": "Factory",
    "title": "Base.:==",
    "category": "method",
    "text": "==(f1::F,f2::F) where F<:AbstractFactory -> Bool\nisequal(f1::F,f2::F) where F<:AbstractFactory -> Bool\n\nOverloaded equivalent operator.\n\n\n\n\n\n"
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
    "location": "man/Utilities/CompositeStructure.html#",
    "page": "Composite structure",
    "title": "Composite structure",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Utilities.CompositeStructure"
},

{
    "location": "man/Utilities/CompositeStructure.html#Composite-structure-1",
    "page": "Composite structure",
    "title": "Composite structure",
    "category": "section",
    "text": "In principle, Julia is not an object-oriented programming language. For example, only abstract types can be inherited so that subtype cannot inherit fields from their parents. Therefore, Julia prefers composition over inheritance. However, to make a new concrete type behaves much alike another one, tedious reputitions of redifining the generic interfaces are usually not avoidable, especially for the basic types in Julia base. In this module, we implement three such composited types, CompositeNTuple, CompositeVector and CompositeDict, for the sake of future usages."
},

{
    "location": "man/Utilities/CompositeStructure.html#CompositeNTuple-1",
    "page": "Composite structure",
    "title": "CompositeNTuple",
    "category": "section",
    "text": "A composite ntuple is a ntuple that is implemented by including an ordinary NTuple as one of its attributes with the name :contents.To take full advantages of the Julia base, the following interfaces are defined:inquiry of info: length, eltype\ncomparison between objects: ==, isequal\nobtainment of old elements: getindex\niteration: iterate, keys, values, pairsComposite ntuples are suited for the situations where other attributes are not affected by the modification of the elements. Note that arithmatic operations and logical operations excluding == and isequal are not supported. Besides, a composite ntuple is not a tuple since Julia has no abstract tuples."
},

{
    "location": "man/Utilities/CompositeStructure.html#CompositeVector-1",
    "page": "Composite structure",
    "title": "CompositeVector",
    "category": "section",
    "text": "A composite vector is a vector that is implemented by including an ordinary Vector as one of its attributes with the name :contents.To take full advantages of the Julia base, the following interfaces are redined:inquiry of info: size, length\ncomparison between objects: ==, isequal\nobtainment of old elements: getindex\nmodification of old elements: setindex!\naddition of new elements: push!, pushfirst!, insert!, append!, prepend!\nremoval of old elements: splice!, deleteat!, pop!, popfirst!, empty!\nconstruction of new objects: empty, reverse, similar\niteration: iterate, keys, values, pairsComposite vectors are suited for the situations where other attributes are not affected by the modification of the elements. Note that arithmatic operations and logical operations excluding == and isequal are not supported."
},

{
    "location": "man/Utilities/CompositeStructure.html#CompositeDict-1",
    "page": "Composite structure",
    "title": "CompositeDict",
    "category": "section",
    "text": "A composite dict is a dict that is implemented by including an ordinary Dict as one of its attributes with the name :contents.To take full advantages of the Julia base, the following interfaces are redined:inquiry of info: isempty, length, haskey, in, hash\ncomparison between objects: ==, isequal\nobtainment of old elements: get, getkey, getindex\nmodification and addition of elements: push!, get!, setindex!\nremoval of old elements: pop!, delete!, empty!\nconstruction of new objects: merge, empty\niteration: iterate, keys, values, pairsAs is similar to composite vectors, composite dicts are suited for the situations where other attributes are not affected by the modification of the elements."
},

{
    "location": "man/Utilities/CompositeStructure.html#Hamiltonian.Utilities.CompositeStructure.CompositeDict",
    "page": "Composite structure",
    "title": "Hamiltonian.Utilities.CompositeStructure.CompositeDict",
    "category": "type",
    "text": "CompositeDict{K,V}\n\nA composite dict is a dict that is implemented by including an ordinary Dict as one of its attributes with the name :contents.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/CompositeStructure.html#Hamiltonian.Utilities.CompositeStructure.CompositeNTuple",
    "page": "Composite structure",
    "title": "Hamiltonian.Utilities.CompositeStructure.CompositeNTuple",
    "category": "type",
    "text": "CompositeNTuple{N,T}\n\nA composite ntuple is a ntuple that is implemented by including an ordinary NTuple as one of its attributes with the name :contents.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/CompositeStructure.html#Hamiltonian.Utilities.CompositeStructure.CompositeVector",
    "page": "Composite structure",
    "title": "Hamiltonian.Utilities.CompositeStructure.CompositeVector",
    "category": "type",
    "text": "CompositeVector{T}\n\nA composite vector is a vector that is implemented by including an ordinary Vector as one of its attributes with the name :contents.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/CompositeStructure.html#Manul-1",
    "page": "Composite structure",
    "title": "Manul",
    "category": "section",
    "text": "Modules=[CompositeStructure]\nOrder=  [:module,:constant,:type,:macro,:function]"
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
    "text": "AbstractTree{N,D} is the abstract type for all concrete trees. By design, it has two type parameters:N: the type of the tree\'s node\nD: the type of the tree\'s dataTo fully utilize the methods designed for a tree structure, in our protocol, a concrete subtype must implement the following methods:inquiry related methods\nroot(tree::AbstractTree{N,D}) where {N,D} -> Union{N,Nothing}\nhaskey(tree::AbstractTree{N,D},node::N) where {N,D} -> Bool\nlength(tree::AbstractTree) -> Int\nparent(tree::AbstractTree{N,D},node::N,superparent::Union{N,Nothing}=nothing) where {N,D} -> Union{N,Nothing}\nchildren(tree::AbstractTree{N,D},node::N) where {N,D} -> Vector{N}\nGet a tree\'s root node (nothing for empty trees)\nGet the number of a tree\'s nodes.\nCheck whether a node is in a tree.\nGet the parent of a tree\'s node or return superparent when the input node is the tree\'s root.\nGet the children of a tree\'s node.\nstructure modification related methods\naddnode!(tree::AbstractTree{N,D},parent::Union{N,Nothing},node::N) where {N,D}\ndeletenode!(tree::AbstractTree{N,D},node::N) where {N,D}\nUpdate the structure of a tree by adding a node. When the parent is nothing, the input tree must be empty and the input node becomes the tree\'s root.\nUpdate the structure of a tree by deleting a node.\nindex related methods\ngetindex(tree::AbstractTree{N,D},node::N) where {N,D} -> D\nsetindex!(tree::AbstractTree{N,D},node::N,data::D) where {N,D}\nGet the data of a tree\'s node\nSet the data of a tree\'s node.Based on these methods, we implement several generic functions for inquiries and manipulationsinquiry for type parameters: keytype, valtype, eltype\nexpansion over nodes/data-records: keys, values, pairs\ninquiry for info of nodes: isleaf, level\ninquiry for nodes: ancestor, descendants, siblings, leaves\nmodification: push!, append!, delete!, empty!And optionally, when a subtype implement the following method,empty(tree::AbstractTree) -> typeof(tree)which constructs an empty tree of the same type with the input one, two more more methods are supported:subtree: Get a subtree starting from a node.\nmove!: Move a subtree to a new position."
},

{
    "location": "man/Utilities/Tree.html#TreeCore-and-SimpleTree-1",
    "page": "Tree",
    "title": "TreeCore and SimpleTree",
    "category": "section",
    "text": "To implement all the prerequisites listed above costs a bit efforts. We provide two lazy ways to get over this:Inheritance AbstractTree with TREECORE::TreeCore as the last attribute\nInclusion an attribute which is an instance of SimpleTree"
},

{
    "location": "man/Utilities/Tree.html#TreeCore-1",
    "page": "Tree",
    "title": "TreeCore",
    "category": "section",
    "text": "TreeCore{N,D}, as the literal meaning indicates, is the core of a tree. It encapsulates all the data structures needed by the default implementation, which constains 4 attributes:root::N: the tree\'s root node\ncontents::Dict{N,D}: the tree\'s (node,data) pairs\nparent::Dict{N,N}: records of the parent of each of the tree\'s nodes\nchildren::Dict{N,Vector{N}}: records of the children of each of the tree\'s nodesAs above, the first lazy way is to include this struct with the special name :TREECORE in your concrete subtype as the last attribute. This process can be even lazier, in that we provide a macro @tree to decorate your \"raw\" struct automatically, e.g.@tree struct SimpleSubTree end\n@tree struct SubTreeWithTreeParameters end {N<:AbstractString,D<:Number}\n@tree struct SubTreeWithCertainTreeParameters end {<:String,<:Int}\n@tree struct SubTreeWithFields info::Vector{Int} end {N<:AbstractString,D<:Number}\n@tree struct SubTreeWithParametricFields{T} info::Vector{T} end {N<:AbstractString,D<:Number}\n@tree struct SubTreeWithOverlappedParametricFields{N} info::Vector{N} end {N<:AbstractString,D<:Number}"
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
    "text": "treedepth\n\nIndicate that the iteration over a tree is depth-first.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Hamiltonian.Utilities.Tree.treewidth",
    "page": "Tree",
    "title": "Hamiltonian.Utilities.Tree.treewidth",
    "category": "constant",
    "text": "treewidth\n\nIndicate that the iteration over a tree is width-first.\n\n\n\n\n\n"
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
    "text": "SimpleTree{N,D}() where {N,D}\n\nThe minimum tree structure that implements all the default tree methods.\n\n\n\n\n\n"
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
    "text": "@tree structdef treeparams::Union{Expr,Nothing}=nothing\n\nDecorate a \"raw\" struct to be a subtype of AbstractTree.\n\nnote: Note\nA \"raw\" struct means:\nIt has no explicit supertype;\nIt has no inner constructor;\nIt has no attribute :TREECORE.\nThe keytype and valtype can be assigned by the argument treeparams in the form {keytype,valtype}.\nWhen the formal argument names of keytype and valtype are not assigned, they can be automatically generated by the functioin gensym. For example, all of the structs after the decration by the following codes\n@tree struct SubTreeWithWrongTypeParameterNames{N} info::Vector{N} end\n@tree struct SubTreeWithWrongTypeParameterNames{N} info::Vector{N} end {::String,::Int}\n@tree struct SubTreeWithWrongTypeParameterNames{N} info::Vector{N} end {<:AbstractString,<:Number}\nwill have three type parameters.\nWhen the formal argument names of keytype and valtype overlap with those of the raw struct type parameters, the duplicates will be considered as the same. For example, the decorated struct SubTreeWithOverlappedParametricFields by the following code\n@tree struct SubTreeWithOverlappedParametricFields{N} info::Vector{N} end {N<:AbstractString,D<:Number}\nonly has two type parameters N<:AbstractString and D<:Number, where the N in the info::Vector{N} is the same N with that in the decorated attribute TREECORE::TreeCore{N,D}.\nWhen the formal argument names of keytype and valtype have no intersection with those of the raw struct type parameters, the type parameters of the decorated struct will be just extended by keytype and valtype. For example, the decorated struct SubTreeWithParametricFields by the following code\n@tree struct SubTreeWithParametricFields{T} info::Vector{T} end {N<:AbstractString,D<:Number}\nhave 3 type parameters, T, N<:AbstractString and D<:Number.\n\n\n\n\n\n"
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
    "text": "==(tc1::TC,tc2::TC) where TC<:TreeCore -> Bool\nisequal(tc1::TC,tc2::TC) where TC<:TreeCore -> Bool\n\nOverloaded equivalent operator.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/Tree.html#Base.:==-Union{Tuple{T}, Tuple{T,T}} where T<:Hamiltonian.Utilities.Tree.AbstractTree",
    "page": "Tree",
    "title": "Base.:==",
    "category": "method",
    "text": "==(t1::T,t2::T) where T<:AbstractTree -> Bool\nisequal(t1::T,t2::T) where T<:AbstractTree -> Bool\n\nOverloaded equivalent operator.\n\n\n\n\n\n"
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
    "location": "man/Utilities/Tree.html#Base.empty-Union{Tuple{AbstractTree{N,D}}, Tuple{D}, Tuple{N}} where D where N",
    "page": "Tree",
    "title": "Base.empty",
    "category": "method",
    "text": "empty(tree::AbstractTree)\n\nConstruct an empty tree of the same type with the input one.\n\n\n\n\n\n"
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
    "location": "man/Utilities/Tree.html#Base.keys-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},TreeIteration}, Tuple{AbstractTree{N,D},TreeIteration,Union{Nothing, N}}} where D where N",
    "page": "Tree",
    "title": "Base.keys",
    "category": "method",
    "text": "keys(tree::AbstractTree{N,D},::TreeDepth,node::Union{N,Nothing}=tree|>root) where {N,D}\nkeys(tree::AbstractTree{N,D},::TreeWidth,node::Union{N,Nothing}=tree|>root) where {N,D}\n\nIterate over a tree\'s nodes starting from a certain node by depth first search or width first search.\n\n\n\n\n\n"
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
    "location": "man/Utilities/Tree.html#Base.pairs-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree,TreeIteration}, Tuple{AbstractTree,TreeIteration,Union{Nothing, N}}} where D where N",
    "page": "Tree",
    "title": "Base.pairs",
    "category": "method",
    "text": "pairs(tree::AbstractTree,::TreeDepth,node::Union{N,Nothing}=tree|>root) where {N,D}\npairs(tree::AbstractTree,::TreeWidth,node::Union{N,Nothing}=tree|>root) where {N,D}\n\nIterate over a tree\'s (node,data) pairs starting from a certain node by depth first search or width first search.\n\n\n\n\n\n"
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
    "location": "man/Utilities/Tree.html#Base.values-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree,TreeIteration}, Tuple{AbstractTree,TreeIteration,Union{Nothing, N}}} where D where N",
    "page": "Tree",
    "title": "Base.values",
    "category": "method",
    "text": "values(tree::AbstractTree,::TreeDepth,node::Union{N,Nothing}=tree|>root) where {N,D}\nvalues(tree::AbstractTree,::TreeWidth,node::Union{N,Nothing}=tree|>root) where {N,D}\n\nIterate over a tree\'s data starting from a certain node by depth first search or width first search.\n\n\n\n\n\n"
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
    "text": "HomoNamedVector is the subtype of [AbstractNamedVector] that of all its fields share the same type. Compared to AbstractNamedVector, one more default method is implemented with HomoNamedVector, i.e. eltype, which returns the type of its fields. This function ensures the type stability of all the methods that involves an iteration of the field values of a named vector. Therefore, homogeneous named vector are usually more efficient than inhomogeneous ones. Use homogeneous ones as much as possible unless the code efficiency does not matter.To subtype HomoNamedVector, all the suggestions mentioned in the previous subsection for AbstractNamedVector also applies. A recommended template for a subtype is[mutable] struct YourNamedVector{T} <: HomoNamedVector{T}\n    filed1::T\n    filed2::T\n    ...\nendWe also provide a macro @homonamedvector to help the definition of concrete homogeneous named vector, where you only need specify the type name, field names, data type and optionally whether the subtype is mutable. For example,@homonamedvector HomoNVWithoutParameter (:scope,:site) Int mutable=true\n@homonamedvector HomoNVWithParameter (:scope,:site) (<:Real) mutable=trueThis macro also integrates the Base.fieldnames function, thus its overloading by hand is on longer needed."
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
    "location": "man/Utilities/NamedVector.html#Base.:<-Tuple{Hamiltonian.Utilities.NamedVector.AbstractNamedVector,Hamiltonian.Utilities.NamedVector.AbstractNamedVector}",
    "page": "Named vector",
    "title": "Base.:<",
    "category": "method",
    "text": "<(nv1::AbstractNamedVector,nv2::AbstractNamedVector) -> Bool\nisless(nv1::AbstractNamedVector,nv2::AbstractNamedVector) -> Bool\n\nCompare two named vectors and judge whether the first is less than the second.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.:==-Tuple{Hamiltonian.Utilities.NamedVector.AbstractNamedVector,Hamiltonian.Utilities.NamedVector.AbstractNamedVector}",
    "page": "Named vector",
    "title": "Base.:==",
    "category": "method",
    "text": "==(nv1::AbstractNamedVector,nv2::AbstractNamedVector) -> Bool\nisequal(nv1::AbstractNamedVector,nv2::AbstractNamedVector) -> Bool\n\nOverloaded equivalent operator. Two named vector are equal to each other if and only if their keys as well as their values are equal to each other.\n\nnote: Note\nIt is not necessary for two named vectors to be of the same concrete type to be equal to each other.\n\n\n\n\n\n"
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
    "text": "keys(nv::AbstractNamedVector) -> NTuple(nv|>fieldcount,Symbol)\n\nIterate over the names.\n\n\n\n\n\n"
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
    "location": "man/Utilities/AlgebraOverField.html#",
    "page": "Algebra over fields",
    "title": "Algebra over fields",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Utilities.AlgebraOverFieldpush!(LOAD_PATH,\"../../../../src/\")\nusing Hamiltonian.Utilities.AlgebraOverField"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Algebra-over-fields-1",
    "page": "Algebra over fields",
    "title": "Algebra over fields",
    "category": "section",
    "text": ""
},

{
    "location": "man/Utilities/AlgebraOverField.html#ID-1",
    "page": "Algebra over fields",
    "title": "ID",
    "category": "section",
    "text": ""
},

{
    "location": "man/Utilities/AlgebraOverField.html#VectorSpace-1",
    "page": "Algebra over fields",
    "title": "VectorSpace",
    "category": "section",
    "text": ""
},

{
    "location": "man/Utilities/AlgebraOverField.html#Element-and-Elements-1",
    "page": "Algebra over fields",
    "title": "Element and Elements",
    "category": "section",
    "text": ""
},

{
    "location": "man/Utilities/AlgebraOverField.html#Hamiltonian.Utilities.AlgebraOverField.VectorSpace",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Utilities.AlgebraOverField.VectorSpace",
    "category": "constant",
    "text": "VectorSpace\n\nThe corresponding vector space of an algebra over a field.\n\nAlias for Union{SimpleVectorSpace,CompositeVectorSpace}.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Hamiltonian.Utilities.AlgebraOverField.VectorSpace-Tuple{Vararg{Hamiltonian.Utilities.AlgebraOverField.SimpleID,N} where N}",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Utilities.AlgebraOverField.VectorSpace",
    "category": "method",
    "text": "VectorSpace(ids::SimpleID...) -> SimpleVectorSpace\nVectorSpace(svses::SimpleVectorSpace...) -> CompositeVectorSpace\n\nGet the corresponding vector space of an algebra over a field.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Hamiltonian.Utilities.AlgebraOverField.CompositeVectorSpace",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Utilities.AlgebraOverField.CompositeVectorSpace",
    "category": "type",
    "text": "CompositeVectorSpace(svses::SimpleVectorSpace...)\n\nThe vector space spanned by the direct product of simple vector spaces.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Hamiltonian.Utilities.AlgebraOverField.Element",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Utilities.AlgebraOverField.Element",
    "category": "type",
    "text": "Element{V<:Number,I<:ID}\n\nAn element of an algebra over a field.\n\nThe first and second attributes of an element must be\n\nvalue::Nuber: the coefficient of the element\nid::ID: the id of the element\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Hamiltonian.Utilities.AlgebraOverField.Elements",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Utilities.AlgebraOverField.Elements",
    "category": "type",
    "text": "Elements{I<:ID,M<:Element} <: AbstractDict{I,M}\n\nAn set of elements of an algebra over a field.\n\nAlias for Dict{I<:ID,M<:Element}. Similar iterms are automatically merged thanks to the id system.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Hamiltonian.Utilities.AlgebraOverField.Elements-Tuple{Any}",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Utilities.AlgebraOverField.Elements",
    "category": "method",
    "text": "Elements(ms)\nElements(ms::Pair{I,M}...) where {I<:ID,M<:Element}\nElements(ms::Element...)\n\nGet the set of elements with similar items merged.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Hamiltonian.Utilities.AlgebraOverField.ID",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Utilities.AlgebraOverField.ID",
    "category": "type",
    "text": "ID(ids::NTuple{N,SimpleID}) where N\nID(ids::SimpleID...)\nID(::Type{SID},attrs::Vararg{NTuple{N},M}) where {SID<:SimpleID,N,M}\n\nThe id system of an algebra over a field.\n\nUsually, a simple id corresponds to a single generator of the algebra while an id corresponds to an element of the algebra.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Hamiltonian.Utilities.AlgebraOverField.SimpleID",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Utilities.AlgebraOverField.SimpleID",
    "category": "type",
    "text": "SimpleID <: AbstractNamedVector\n\nA simple id is the building block of the id system of an algebra over a field.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Hamiltonian.Utilities.AlgebraOverField.SimpleVectorSpace",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Utilities.AlgebraOverField.SimpleVectorSpace",
    "category": "type",
    "text": "SimpleVectorSpace(ids::SimpleID...)\n\nThe vector space spanned by a set of bases specified by their ids.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Hamiltonian.Utilities.:⊕-Union{Tuple{I}, Tuple{I,I}} where I<:Hamiltonian.Utilities.AlgebraOverField.SimpleID",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Utilities.:⊕",
    "category": "method",
    "text": "⊕(id1::I,id2::I) where {I<:SimpleID} -> SimpleVectorSpace{I}\n⊕(id::I,svs::SimpleVectorSpace{I}) where {I<:SimpleID} -> SimpleVectorSpace{I}\n⊕(svs::SimpleVectorSpace{I},id::I) where {I<:SimpleID} -> SimpleVectorSpace{I}\n⊕(svs1::SVS,svs2::SVS) where {SVS<:SimpleVectorSpace} -> SVS\n\nGet the direct sum of bases or simple vector spaces.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Hamiltonian.Utilities.:⊗-Tuple{Hamiltonian.Utilities.AlgebraOverField.SimpleID,Hamiltonian.Utilities.AlgebraOverField.SimpleID}",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Utilities.:⊗",
    "category": "method",
    "text": "⊗(sid1::SimpleID,sid2::SimpleID) -> ID\n⊗(sid::SimpleID,cid::ID) -> ID\n⊗(cid::ID,sid::SimpleID) -> ID\n⊗(cid1::ID,cid2::ID) -> ID\n\nGet the direct product of the id system.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Hamiltonian.Utilities.:⊗-Tuple{Hamiltonian.Utilities.AlgebraOverField.SimpleVectorSpace,Hamiltonian.Utilities.AlgebraOverField.SimpleVectorSpace}",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Utilities.:⊗",
    "category": "method",
    "text": "⊗(svs1::SimpleVectorSpace,svs2::SimpleVectorSpace) -> CompositeVectorSpace\n⊗(svs::SimpleVectorSpace,cvs::CompositeVectorSpace) -> CompositeVectorSpace\n⊗(cvs::CompositeVectorSpace,svs::SimpleVectorSpace) -> CompositeVectorSpace\n⊗(cvs1::CompositeVectorSpace,cvs2::CompositeVectorSpace) -> CompositeVectorSpace\n\nGet the direct product of simple vector spaces or composite vector spaces.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Hamiltonian.Utilities.AlgebraOverField.add!-Tuple{Dict{I,M} where M<:Hamiltonian.Utilities.AlgebraOverField.Element where I<:Hamiltonian.Utilities.AlgebraOverField.ID}",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Utilities.AlgebraOverField.add!",
    "category": "method",
    "text": "add!(ms::Elements) -> typeof(ms)\nadd!(ms::Elements,m::Element) -> typeof(ms)\nadd!(ms::Elements,mms::Elements) -> typeof(ms)\n\nGet the inplace addition of elements to a set.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Hamiltonian.Utilities.AlgebraOverField.idtype-Union{Tuple{Type{#s68} where #s68<:(Element{V,I,N} where N)}, Tuple{I}, Tuple{V}} where I where V",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Utilities.AlgebraOverField.idtype",
    "category": "method",
    "text": "idtype(::Type{<:Element{V,I}}) where {V,I}\nidtype(m::Element)\n\nThe type of the id of an element.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Hamiltonian.Utilities.AlgebraOverField.sub!-Tuple{Dict{I,M} where M<:Hamiltonian.Utilities.AlgebraOverField.Element where I<:Hamiltonian.Utilities.AlgebraOverField.ID}",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Utilities.AlgebraOverField.sub!",
    "category": "method",
    "text": "sub!(ms::Elements) -> typeof(ms) -> typeof(ms)\nsub!(ms::Elements,m::Element) -> typeof(ms)\nsub!(ms::Elements,mms::Elements) -> typeof(ms)\n\nGet the inplace subtraction of elements from a set.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Hamiltonian.Utilities.dimension-Tuple{Hamiltonian.Utilities.AlgebraOverField.SimpleVectorSpace}",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Utilities.dimension",
    "category": "method",
    "text": "dimension(svs::SimpleVectorSpace) -> Int\ndimension(cvs::CompositeVectorSpace) -> Int\n\nGet the dimension of a vector space.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Hamiltonian.Utilities.rank-Union{Tuple{Type{#s68} where #s68<:Element{V,I,N}}, Tuple{N}, Tuple{I}, Tuple{V}} where N where I where V",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Utilities.rank",
    "category": "method",
    "text": "rank(::Type{<:Element{V,I,N}}) where {V,I,N} -> Int\nrank(m::Element) -> Int\n\nGet the rank of an element.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Hamiltonian.Utilities.rank-Union{Tuple{Type{#s68} where #s68<:ID{N,I}}, Tuple{I}, Tuple{N}} where I where N",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Utilities.rank",
    "category": "method",
    "text": "rank(::Type{<:ID{N,I}}) where {N,I} -> Int\nrank(id::ID) -> Int\n\nGet the rank of a composite id.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Base.:*-Tuple{Number,Hamiltonian.Utilities.AlgebraOverField.Element}",
    "page": "Algebra over fields",
    "title": "Base.:*",
    "category": "method",
    "text": "*(factor::Number,m::Element) -> Element\n*(m::Element,factor::Number) -> Element\n*(m1::Element,m2::Element) -> Element\n*(factor::Number,ms::Elements) -> Elements\n*(ms::Elements,factor::Number) -> Elements\n*(m::Element,ms::Elements) -> Elements\n*(ms::Elements,m::Element) -> Elements\n*(ms1::Elements,ms2::Elements) -> Elements\n\nOverloaded * operator for element-scalar multiplications and element-element multiplications of an algebra over a field.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Base.:+-Tuple{Hamiltonian.Utilities.AlgebraOverField.Element}",
    "page": "Algebra over fields",
    "title": "Base.:+",
    "category": "method",
    "text": "+(m::Element) -> typeof(m)\n+(ms::Elements) -> typeof(ms)\n+(ms::Elements,m::Element) -> Elements\n+(m1::Element,m2::Element) -> Elements\n+(m::Element,ms::Elements) -> Elements\n+(ms1::Elements,ms2::Elements) -> Elements\n\nOverloaded + operator between elements of an algebra over a field.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Base.:--Tuple{Hamiltonian.Utilities.AlgebraOverField.Element}",
    "page": "Algebra over fields",
    "title": "Base.:-",
    "category": "method",
    "text": "-(m::Element) -> typeof(m)\n-(ms::Elements) -> typeof(ms)\n-(m1::Element,m2::Element) -> Elements\n-(m::Element,ms::Elements) -> Elements\n-(ms::Elements,m::Element) -> Elements\n-(ms1::Elements,ms2::Elements) -> Elements\n\nOverloaded - operator between elements of an algebra over a field.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Base.:/-Tuple{Hamiltonian.Utilities.AlgebraOverField.Element,Number}",
    "page": "Algebra over fields",
    "title": "Base.:/",
    "category": "method",
    "text": "/(m::Element,factor::Number)\n/(ms::Elements,factor::Number)\n\nOverloaded / operator for element-sclar division of an algebra over a field.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Base.:==-Union{Tuple{M}, Tuple{M,M}} where M<:Hamiltonian.Utilities.AlgebraOverField.Element",
    "page": "Algebra over fields",
    "title": "Base.:==",
    "category": "method",
    "text": "==(m1::M,m2::M) where M<:Element -> Bool\nisequal(m1::M,m2::M) where M<:Element -> Bool\n\nCompare two elements and judge whether they are equal to each other.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Base.getproperty-Tuple{Hamiltonian.Utilities.AlgebraOverField.ID,Symbol}",
    "page": "Algebra over fields",
    "title": "Base.getproperty",
    "category": "method",
    "text": "getproperty(cid::ID,name::Symbol)\n\nGet the property of a composite id.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Base.hash-Tuple{Hamiltonian.Utilities.AlgebraOverField.ID,UInt64}",
    "page": "Algebra over fields",
    "title": "Base.hash",
    "category": "method",
    "text": "hash(cid::ID,h::UInt)\n\nHash a composite id.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Base.isless-Tuple{Hamiltonian.Utilities.AlgebraOverField.ID,Hamiltonian.Utilities.AlgebraOverField.ID}",
    "page": "Algebra over fields",
    "title": "Base.isless",
    "category": "method",
    "text": "isless(cid1::ID,cid2::ID) -> Bool\n<(cid1::ID,cid2::ID) -> Bool\n\nCompare two ids and judge whether the first is less than the second.\n\nWe assume that ids with smaller ranks are always less than those with higher ranks. If two ids are of the same rank, the comparison goes just like that between tuples.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Base.propertynames-Union{Tuple{Type{I}}, Tuple{I}, Tuple{Type{I},Bool}} where I<:Hamiltonian.Utilities.AlgebraOverField.ID",
    "page": "Algebra over fields",
    "title": "Base.propertynames",
    "category": "method",
    "text": "propertynames(::Type{I},private::Bool=false) where I<:ID -> Tuple\n\nGet the property names of a composite id.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Base.replace-Tuple{Hamiltonian.Utilities.AlgebraOverField.Element}",
    "page": "Algebra over fields",
    "title": "Base.replace",
    "category": "method",
    "text": "replace(m::Element;kwargs...) -> typeof(m)\n\nReturn a copy of a concrete Element with some of the field values replaced by the keyword arguments.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Base.show-Tuple{IO,Hamiltonian.Utilities.AlgebraOverField.ID}",
    "page": "Algebra over fields",
    "title": "Base.show",
    "category": "method",
    "text": "show(io::IO,cid::ID)\n\nShow a composite id.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Base.valtype-Union{Tuple{Type{#s68} where #s68<:(Element{V,I,N} where N where I<:ID)}, Tuple{V}} where V",
    "page": "Algebra over fields",
    "title": "Base.valtype",
    "category": "method",
    "text": "valtype(::Type{<:Element{V}}) where {V}\nvaltype(m::Element)\n\nGet the type of the value of an element.\n\nThe result is also the type of the field over which the algebra is defined.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Base.zero-Tuple{Dict{I,M} where M<:Hamiltonian.Utilities.AlgebraOverField.Element where I<:Hamiltonian.Utilities.AlgebraOverField.ID}",
    "page": "Algebra over fields",
    "title": "Base.zero",
    "category": "method",
    "text": "zero(ms::Elements) -> typeof(ms)\nzero(::Type{Elements{I,M}}) where {I,M} -> Elements{I,M}\n\nGet a zero set of elements.\n\nA zero set of elements is defined to be the empty one.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/AlgebraOverField.html#Manual-1",
    "page": "Algebra over fields",
    "title": "Manual",
    "category": "section",
    "text": "Modules=[AlgebraOverField]\nOrder=  [:module,:constant,:type,:macro,:function]"
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
    "text": "qnsbruteforce\n\nIndicate that decompose uses the brute force method.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.qnscompression",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.qnscompression",
    "category": "constant",
    "text": "qnscompression\n\nIndicate that findall and permute use the compressed contents.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.qnscontents",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.qnscontents",
    "category": "constant",
    "text": "qnscontents\n\nIndicate that expand uses the compressed/expanded contents.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.qnscounts",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.qnscounts",
    "category": "constant",
    "text": "qnscounts\n\nIndicate that methods with QuantumNumbers use the count number of the compressed contents.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.qnsexpansion",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.qnsexpansion",
    "category": "constant",
    "text": "qnsexpansion\n\nIndicate that findall and permute use the expanded contents.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.qnsindices",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.qnsindices",
    "category": "constant",
    "text": "qnsindices\n\nIndicate that expand uses the indices of the compressed/expanded contents.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.qnsindptr",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.qnsindptr",
    "category": "constant",
    "text": "qnsindptr\n\nIndicate that methods with QuantumNumbers use the index pointer of the compressed contents.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.qnsmontecarlo",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.qnsmontecarlo",
    "category": "constant",
    "text": "qnsmontecarlo\n\nIndicate that decompose uses the Monte Carlo method.\n\n\n\n\n\n"
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
    "text": "QuantumNumbers(form::Char,contents::Vector{<:AbstractQuantumNumber},counts::Vector{Int},::QnsCounts)\nQuantumNumbers(form::Char,contents::Vector{<:AbstractQuantumNumber},indptr::Vector{Int},::QnsIndptr)\n\nThe whole quantum numbers of the total bases of a Hilbert space.\n\nThe default constructors construct a QuantumNumbers from a vector of concrete quantum numbers and an vector containing their counts or indptr.\n\n\n\n\n\n"
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
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.:⊕-Tuple{Vararg{Hamiltonian.Utilities.QuantumNumber.AbstractQuantumNumber,N} where N}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.:⊕",
    "category": "method",
    "text": "⊕(qns::AbstractQuantumNumber...) -> QuantumNumbers{qns|>eltype}\n⊕(qnses::QuantumNumbers...) -> qnses|>eltype\n\nGet the direct sum of some AbstractQuantumNumbers or QuantumNumberses.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.:⊗-Tuple{Vararg{Hamiltonian.Utilities.QuantumNumber.AbstractQuantumNumber,N} where N}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.:⊗",
    "category": "method",
    "text": "⊗(qns::AbstractQuantumNumber...) -> eltype(qns)\n⊗(qnses::QuantumNumbers...) -> eltype(qnses)\n\nGet the direct product of some AbstractQuantumNumbers or QuantumNumberses.\n\n\n\n\n\n"
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
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.decompose-Union{Tuple{QN}, Tuple{N}, Tuple{Tuple{Vararg{QuantumNumbers{QN},N}},QN,Tuple{Vararg{Int64,N}},QnsBruteForce}} where QN<:Hamiltonian.Utilities.QuantumNumber.AbstractQuantumNumber where N",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.decompose",
    "category": "method",
    "text": "decompose(qnses::NTuple{N,QuantumNumbers{QN}},target::QN,signs::NTuple{N,Int},::QnsBruteForce;nmax::Int=20) where {N,QN<:AbstractQuantumNumber} -> Vector{NTuple{N,Int}}\ndecompose(qnses::NTuple{N,QuantumNumbers{QN}},target::QN,signs::NTuple{N,Int},::QnsMonteCarlo;nmax::Int=20) where {N,QN<:AbstractQuantumNumber} -> Vector{NTuple{N,Int}}\n\nFind a couple of decompositions of target with respect to qnses.\n\nnote: Note\nA tuple of integers (i₁,i₂,...) is called a decomposition of a given target with respect to the given qnses if and only if they satisfy the \"decomposition rule\":sum_textj textsignstextjtimestextqnsestextjtexti_textj==texttargetThis equation is in fact a kind of a set of restricted linear Diophantine equations. Indeed, our quantum numbers are always discrete Abelian ones and all instances of a concrete AbstractQuantumNumber forms a module over the ring of integers. Therefore, each quantum number can be represented as a integral multiple of the unit element of the Abelian module, which results in the final reduction of the above equation to a set of linear Diophantine equations. Then finding a decomposition is equivalent to find a solution of the reduced linear Diophantine equations, with the restriction that the quantum numbers constructed from the solution should be in the corresponding qnses. Here we provide two methods to find such decompositions, one is by brute force (qnsbruteforce case), and the other is by Monte Carlo simultatioins (qnsmontecarlo case).\n\n\n\n\n\n"
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
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.QuantumNumber.toordereddict-Tuple{Hamiltonian.Utilities.QuantumNumber.QuantumNumbers,Hamiltonian.Utilities.QuantumNumber.QnsIndptr}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.QuantumNumber.toordereddict",
    "category": "method",
    "text": "toordereddict(qns::QuantumNumbers,::QnsIndptr) -> OrderedDict{qns|>eltype,UnitRange{Int}}\ntoordereddict(qns::QuantumNumbers,::QnsCounts) -> OrderedDict{qns|>eltype,Int}\n\nConvert a QuantumNumbers to an ordered dict.\n\n\n\n\n\n"
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
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.expand-Tuple{Hamiltonian.Utilities.QuantumNumber.QuantumNumbers,Hamiltonian.Utilities.QuantumNumber.QnsContents}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.expand",
    "category": "method",
    "text": "expand(qns::QuantumNumbers,::QnsContents) -> Vector{qns|>eltype}\nexpand(qns::QuantumNumbers,::QnsIndices) -> Vector{Int}\n\nExpand the contents (qnscontents case) or indices (qnsindices case) of a QuantumNumbers to the uncompressed form.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Hamiltonian.Utilities.permute-Tuple{Hamiltonian.Utilities.QuantumNumber.QuantumNumbers,Array{Int64,1},Hamiltonian.Utilities.QuantumNumber.QnsCompression}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Utilities.permute",
    "category": "method",
    "text": "permute(qns::QuantumNumbers,permutation::Vector{Int},::QnsCompression) -> QuantumNumbers\npermute(qns::QuantumNumbers,permutation::Vector{Int},::QnsExpansion) -> QuantumNumbers\n\nReorder the quantum numbers contained in a QuantumNumbers with a permutation and return the new one.\n\nFor qnscompression case, the permutation is for the compressed contents of the original QuantumNumbers while for qnsexpansion case, the permutation is for the expanded contents of the original QuantumNumbers.\n\n\n\n\n\n"
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
    "text": "+(qn::AbstractQuantumNumber) -> typeof(qn)\n+(qn::QN,qns::QN...) where QN<:AbstractQuantumNumber -> QN\n+(qns::QuantumNumbers) -> QuantumNumbers\n+(qn::QN,qns::QuantumNumbers{QN}) where QN<:AbstractQuantumNumber -> QuantumNumbers{QN}\n+(qns::QuantumNumbers{QN},qn::QN) where QN<:AbstractQuantumNumber -> QuantumNumbers{QN}\n\nOverloaded + operator for AbstractQuantumNumber and QuantumNumbers.\n\nnote: Note\nThe addition between a QuantumNumbers and an AbstractQuantumNumber is just a global shift of the contents of the QuantumNumbers by the AbstractQuantumNumber, therefore, the result is a QuantumNumbers.\n+ cannot be used between two QuantumNumbers because the result is ambiguous. Instead, use ⊕ for direct sum and ⊗ for direct product.\nTo ensure type stability, two AbstractQuantumNumber can be added together if and only if they are of the same type.\nSimilarly, a AbstractQuantumNumber and a QuantumNumbers can be added together if and only if the former\'s type is the same with the latter\'s eltype.\n\n\n\n\n\n"
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
    "text": "==(qns1::QuantumNumbers,qns2::QuantumNumbers) -> Bool\nisequal(qns1::QuantumNumbers,qns2::QuantumNumbers) -> Bool\n\nOverloaded equivalent operator. Two QuantumNumberses are equal to each other if and only if both their contentses and indptrs are elementwise equal to each other.\n\nnote: Note\nIt is not necessary for two QuantumNumberses to have the same eltype nor the same form to be equal to each other.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Base.:^-Tuple{Hamiltonian.Utilities.QuantumNumber.AbstractQuantumNumber,Integer}",
    "page": "Quantum numbers",
    "title": "Base.:^",
    "category": "method",
    "text": "^(qn::AbstractQuantumNumber,factor::Integer) -> typeof(qn)\n^(qns::QuantumNumbers,factor::Integer) -> QuantumNumbers\n\nOverloaded ^ operator for AbstractQuantumNumber and QuantumNumbers. This operation translates into the direct product of factor copies of qn or qns.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Base.eltype-Union{Tuple{Type{#s68} where #s68<:QuantumNumbers{QN}}, Tuple{QN}} where QN",
    "page": "Quantum numbers",
    "title": "Base.eltype",
    "category": "method",
    "text": "eltype(::Type{<:QuantumNumbers{QN}}) where QN\neltype(qns::QuantumNumbers)\n\nGet the type of the concrete AbstractQuantumNumber contained in a QuantumNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Base.filter-Union{Tuple{QN}, Tuple{QN,QuantumNumbers{QN}}} where QN<:Hamiltonian.Utilities.QuantumNumber.AbstractQuantumNumber",
    "page": "Quantum numbers",
    "title": "Base.filter",
    "category": "method",
    "text": "filter(target::QN,qns::QuantumNumbers{QN}) where QN<:AbstractQuantumNumber -> QuantumNumbers{QN}\nfilter(targets::NTuple{N,QN},qns::QuantumNumbers{QN}) where {N,QN<:AbstractQuantumNumber} -> QuantumNumbers{QN}\n\nFind a subset of a QuantumNumbers by picking out the quantum numbers in targets.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Base.findall-Union{Tuple{QN}, Tuple{QN,QuantumNumbers{QN},QnsCompression}} where QN<:Hamiltonian.Utilities.QuantumNumber.AbstractQuantumNumber",
    "page": "Quantum numbers",
    "title": "Base.findall",
    "category": "method",
    "text": "findall(target::QN,qns::QuantumNumbers{QN},::QnsCompression) where QN<:AbstractQuantumNumber -> Vector{Int}\nfindall(target::QN,qns::QuantumNumbers{QN},::QnsExpansion) where QN<:AbstractQuantumNumber -> Vector{Int}\nfindall(targets::NTuple{N,QN},qns::QuantumNumbers{QN},::QnsCompression) where {N,QN<:AbstractQuantumNumber} -> Vector{Int}\nfindall(targets::NTuple{N,QN},qns::QuantumNumbers{QN},::QnsExpansion) where {N,QN<:AbstractQuantumNumber} -> Vector{Int}\n\nFind all the indices of the target quantum numbers in the contents (qnscompression case) or the expansion (qnsexpansion case) of a QuantumNumbers.\n\n\n\n\n\n"
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
    "text": "kron(::Type{QN},qn1::AbstractQuantumNumber,qn2::AbstractQuantumNumber) where QN<:AbstractQuantumNumber -> QN\nkron(qns::Vararg{<:AbstractQuantumNumber,N};signs::NTuple{N,Int}=ntuple(i->1,N)) where N -> eltype(qns)\nkron(qnses::Vararg{QuantumNumbers{QN},N};signs::NTuple{N,Int}=ntuple(i->1,N)) where {N,QN<:AbstractQuantumNumber} -> QuantumNumbers{QN}\n\nGet the direct product of some AbstractQuantumNumbers or QuantumNumberses.\n\nnote: Note\nPhysically, the direct product of a couple of AbstractQuantumNumbers or QuantumNumberses are defined by the direct product of the bases of the Hilbert spaces they represent. Therefore, QuantumNumbers with differenct types or QuantumNumberses with differenct eltypes are allowed to be direct producted in principle. However, for simplicity, we only implement a method which handle the situation of two AbstractQuantumNumbers with differenct types. The type of the result should be provided as the first parameter. Note that in this situation, the fieldnames and periods of the result type must be exactly equal to the flattened fieldnames and periods of the two input AbstractQuantumNumbers, which means, even the order of the input AbstractQuantumNumbers matters.\nApparently, the dimension of the result equals the product of those of the inputs. Therefore, the direct product of AbstractQuantumNumbers is also a AbstractQuantumNumber since its dimension is still one.\nFor other situations except the one mentioned in Note.1, the input AbstractQuantumNumbers or QuantumNumberses must be homogenous. Meanwhile, signs can also be provided for these situations. Note that each quantum number in the contents of the result is obtained by a summation of the corresponding quanum numbers out of the inputs with the correct signs. This is a direct observation of the Abelian nature of our quantum numbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Base.length-Tuple{Hamiltonian.Utilities.QuantumNumber.QuantumNumbers}",
    "page": "Quantum numbers",
    "title": "Base.length",
    "category": "method",
    "text": "length(qns::QuantumNumbers) -> Int\n\nGet the number of unduplicate qunatum numbers in the QuantumNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Base.pairs-Tuple{Hamiltonian.Utilities.QuantumNumber.QuantumNumbers,Union{QnsCounts, QnsIndptr}}",
    "page": "Quantum numbers",
    "title": "Base.pairs",
    "category": "method",
    "text": "pairs(qns::QuantumNumbers,choice::Union{QnsIndptr,QnsCounts})\n\nIterate over the AbstractQuantumNumber=>slice or AbstractQuantumNumber=>count pairs.\n\n\n\n\n\n"
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
    "location": "man/Utilities/QuantumNumber.html#Base.union-Union{Tuple{Vararg{AbstractQuantumNumber,N}}, Tuple{N}} where N",
    "page": "Quantum numbers",
    "title": "Base.union",
    "category": "method",
    "text": "union(qns::Vararg{<:AbstractQuantumNumber,N};signs::NTuple{N,Int}=ntuple(i->1,N)) where N -> QuantumNumbers\nunion(qnses::Vararg{QuantumNumbers{QN},N};signs::NTuple{N,Int}=ntuple(i->1,N)) where {N,QN<:AbstractQuantumNumber} -> QuantumNumbers{QN}\n\nGet the direct sum of some AbstractQuantumNumbers or QuantumNumberses.\n\nnote: Note\nPhysically, the direct sum of a couple of AbstractQuantumNumbers or QuantumNumberses is defined by the direct sum of the bases of the Hilbert spaces they represent. Therefore, the input AbstractQuantumNumbers or QuantumNumberses must be homogenous. Inhomogenous \'AbstractQuantumNumber\'s must be direct producted first to ensure homogenity before the direct sum.\nApparently, the dimension of the result equals the summation of those of the inputs, which means, even for AbstractQuantumNumbers, the result will be naturally a QuantumNumbers because the dimension of the result is larger than 1.\nSigns of AbstractQuantumNumbers or QuantumNumberses can be provided when getting their direct sums.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#Base.values-Tuple{Hamiltonian.Utilities.QuantumNumber.QuantumNumbers,Hamiltonian.Utilities.QuantumNumber.QnsIndptr}",
    "page": "Quantum numbers",
    "title": "Base.values",
    "category": "method",
    "text": "values(qns::QuantumNumbers,::QnsIndptr)\nvalues(qns::QuantumNumbers,::QnsCounts)\n\nIterate over the slices/counts of the QuantumNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/QuantumNumber.html#qnmanual-1",
    "page": "Quantum numbers",
    "title": "Manual",
    "category": "section",
    "text": "Modules=[QuantumNumber]\nOrder=  [:module,:constant,:type,:macro,:function]"
},

{
    "location": "man/Essentials/Introduction.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Essentials"
},

{
    "location": "man/Essentials/Introduction.html#Introduction-1",
    "page": "Introduction",
    "title": "Introduction",
    "category": "section",
    "text": "Essentials of the Hamiltonian package, which defines all the imported constants, types and functions when using import Hamiltonian or using Hamiltonian. Note that this submodule depends on the Utilities submodule although the variables in the latter are not exported to the scope of Hamiltonian by default.Pages=  [\n        \"Spatial.md\",\n        \"DegreeOfFreedom.md\",\n        \"FockPackage.md\",\n        ]\nDepth=2"
},

{
    "location": "man/Essentials/Spatial.html#",
    "page": "Spatial",
    "title": "Spatial",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Essentials.Spatialpush!(LOAD_PATH,\"../../../../src/\")\nusing Hamiltonian"
},

{
    "location": "man/Essentials/Spatial.html#Spatial-1",
    "page": "Spatial",
    "title": "Spatial",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/Spatial.html#AbstractBond-1",
    "page": "Spatial",
    "title": "AbstractBond",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/Spatial.html#Point-1",
    "page": "Spatial",
    "title": "Point",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/Spatial.html#Bond-1",
    "page": "Spatial",
    "title": "Bond",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/Spatial.html#AbstractLattice-1",
    "page": "Spatial",
    "title": "AbstractLattice",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/Spatial.html#Lattice-1",
    "page": "Spatial",
    "title": "Lattice",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/Spatial.html#SuperLattice-1",
    "page": "Spatial",
    "title": "SuperLattice",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/Spatial.html#Cylinder-1",
    "page": "Spatial",
    "title": "Cylinder",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.acrossbonds",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.acrossbonds",
    "category": "constant",
    "text": "acrossbonds\n\nIndicate that bonds across the unitcell are inquired, which are in fact those across the periodic boundaries.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.insidebonds",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.insidebonds",
    "category": "constant",
    "text": "insidebonds\n\nIndicate that bonds inside the unitcell are inquired, which do not contain those across the periodic boundaries.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.interbonds",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.interbonds",
    "category": "constant",
    "text": "interbonds\n\nIndicate that bonds inter the sublattices are inquired.\n\nnotes: Notes\nThese bonds do not contain those accorss the periodic boundaries.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.intrabonds",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.intrabonds",
    "category": "constant",
    "text": "intrabonds\n\nIndicate that bonds intra the sublattices are inquired.\n\nnotes: Notes\nThese bonds do not contain those accorss the periodic boundaries.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.zerothbonds",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.zerothbonds",
    "category": "constant",
    "text": "zerothbonds\n\nIndicate that zeroth bonds, i.e. the points are inquired.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.AbstractBond",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.AbstractBond",
    "category": "type",
    "text": "AbstractBond{R,P<:PID,N}\n\nAbstract bond.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.AbstractLattice",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.AbstractLattice",
    "category": "type",
    "text": "AbstractLattice{P<:PID,N}\n\nAbstract type for all lattices.\n\nIt should have the following attributes\n\nname::String: the name of the lattice\npids::Vector{P}: the pids of the lattice\nrcoords::Matrix{Float}: the rcoords of the lattice\nicoords::Matrix{Float}: the icoords of the lattice\nvectors::Vector{SVector{N,Float}}: the translation vectors of the lattice\nreciprocals::Vector{SVector{N,Float}}: the reciprocals of the lattice\nneighbors::Dict{Int,Float}: the order-distance map of the nearest neighbors of the lattice\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.AbstractLatticeIndex",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.AbstractLatticeIndex",
    "category": "type",
    "text": "AbstractLatticeIndex{I<:Union{<:PID,Int}}\n\nAbstract index type for a lattice.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.Bond",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.Bond",
    "category": "type",
    "text": "Bond(neighbor::Int,spoint::Point,epoint::Point)\n\nA bond in a lattice.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.Cylinder",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.Cylinder",
    "category": "type",
    "text": "Cylinder{P}(    name::String,\n                block::AbstractMatrix{<:Real},\n                translation::AbstractVector{<:Real};\n                vector::Union{AbstractVector{<:Real},Nothing}=nothing,\n                neighbors::Union{Dict{Int,<:Real},Int}=1,\n                ) where P<:PID\n\nCylinder of 1d and quasi 2d lattices.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.Cylinder-Tuple",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.Cylinder",
    "category": "method",
    "text": "(cylinder::Cylinder)(scopes::Any...;coordination::Int=8) -> Lattice\n\nConstruct a lattice from a cylinder with the assigned scopes.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.ICoordIndex",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.ICoordIndex",
    "category": "type",
    "text": "ICoordIndex(index::Union{<:PID,Int})\n\nIndex for getting an icoord of a lattice.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.Lattice",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.Lattice",
    "category": "type",
    "text": "Lattice(    name::String,\n            pids::Vector{<:PID},\n            rcoords::AbstractMatrix{<:Real};\n            icoords::AbstractMatrix{<:Real}=SMatrix{0,0,Float}(),\n            vectors::AbstractVector{<:AbstractVector{<:Real}}=SVector{0,SVector{size(rcoords,1),Float}}(),\n            neighbors::Union{Dict{Int,<:Real},Int}=1,\n            coordination::Int=8\n        )\nLattice(    name::String,\n            points::AbstractVector{<:Point};\n            vectors::AbstractVector{<:AbstractVector{<:Real}}=SVector{0,SVector{points|>eltype|>dimension,Float}}(),\n            neighbors::Union{Dict{Int,<:Real},Int}=1,\n            coordination::Int=8\n            )\nLattice(    name::String,\n            sublattices::AbstractVector{<:AbstractLattice};\n            vectors::AbstractVector{<:AbstractVector{<:Real}}=SVector{0,SVector{sublattices|>eltype|>dimension,Float}}(),\n            neighbors::Union{Dict{Int,<:Real},Int}=1,\n            coordination::Int=8\n            )\n\nSimplest lattice.\n\nA simplest lattice can be construted from its contents, i.e. pids, rcoords and icoords, or from a couple of points, or from a couple of sublattices.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.Link",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.Link",
    "category": "type",
    "text": "Link(neighbor::Int,sindex::Int,eindex::Int,disp::AbstractVector{<:Real})\n\nA link in a lattice.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.PID",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.PID",
    "category": "type",
    "text": "PID(scope,site::Int)\nPID(;scope=\"tz\",site::Int=1)\n\nThe id of a point.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.Point",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.Point",
    "category": "type",
    "text": "Point(pid::PID,rcoord::SVector{N,<:Real},icoord::SVector{N,<:Real}) where N\nPoint(pid::PID,rcoord::NTuple{N,<:Real},icoord::NTuple{N,<:Real}=ntuple(i->0.0,N)) where N\nPoint(pid::PID,rcoord::AbstractVector{<:Real},icoord::AbstractVector{<:Real}=zero(SVector{length(rcoord),Float}))\n\nLabeled point.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.PointIndex",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.PointIndex",
    "category": "type",
    "text": "PointIndex(index::Union{<:PID,Int})\n\nIndex for getting a point of a lattice.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.RCoordIndex",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.RCoordIndex",
    "category": "type",
    "text": "RCoordIndex(index::Union{<:PID,Int})\n\nIndex for getting a rcoord of a lattice.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.SuperLattice",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.SuperLattice",
    "category": "type",
    "text": "SuperLattice(   name::String,\n                sublattices::AbstractVector{<:AbstractLattice};\n                vectors::AbstractVector{<:AbstractVector{<:Real}}=SVector{0,SVector{sublattices|>eltype|>dimension,Float}}(),\n                neighbors::Dict{Int,<:Real}=Dict{Int,Float}()\n                )\n\nSuperLattice that is composed of serveral sublattices.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.azimuth-Tuple{AbstractArray{#s199,1} where #s199<:Real}",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.azimuth",
    "category": "method",
    "text": "azimuth(v::AbstractVector{<:Real}) -> Float\n\nGet the azimuth angle in radians of a vector.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.azimuthd-Tuple{AbstractArray{#s199,1} where #s199<:Real}",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.azimuthd",
    "category": "method",
    "text": "azimuthd(v::AbstractVector{<:Real}) -> Float\n\nGet the azimuth angle in degrees of a vector.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.bonds-Tuple{AbstractLattice}",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.bonds",
    "category": "method",
    "text": "bonds(lattice::AbstractLattice) -> Vector{AbstractBond}\nbonds(lattice::AbstractLattice,::ZerothBonds) -> Vector{Point}\nbonds(lattice::AbstractLattice,::InsideBonds) -> Vector{Bond}\nbonds(lattice::AbstractLattice,::AcrossBonds) -> Vector{Bond}\n\nGet the bonds of a lattice.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.bonds-Tuple{SuperLattice,Hamiltonian.Essentials.Spatial.IntraBonds}",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.bonds",
    "category": "method",
    "text": "bonds(lattice::SuperLattice,::IntraBonds) -> Vector{Bond}\nbonds(lattice::SuperLattice,::InterBonds) -> Vector{Bond}\n\nGet the bonds of a superlattice.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.distance-Tuple{AbstractArray{#s188,1} where #s188<:Real,AbstractArray{#s187,1} where #s187<:Real}",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.distance",
    "category": "method",
    "text": "distance(p1::AbstractVector{<:Real},p2::AbstractVector{<:Real}) -> Float\n\nGet the distance between two points.\n\nnote: Note\nCompared to norm(p1-p2), this function avoids the memory allocation for p1-p2, thus is more efficient.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.icoord-Tuple{Bond}",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.icoord",
    "category": "method",
    "text": "icoord(bond::Bond) -> SVector\n\nGet the icoord of the bond.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.interlinks-Tuple{AbstractArray{#s209,2} where #s209<:Real,AbstractArray{#s208,2} where #s208<:Real,Dict{Int64,Float64}}",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.interlinks",
    "category": "method",
    "text": "interlinks(cluster1::AbstractMatrix{<:Real},cluster2::AbstractMatrix{<:Real},neighbors::Dict{Int,Float}) -> Vector{Link}\n\nUse kdtree to get the intercluster nearest neighbors.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.intralinks-Union{Tuple{N}, Tuple{AbstractArray{#s204,2} where #s204<:Real,AbstractArray{#s203,1} where #s203<:(AbstractArray{#s202,1} where #s202<:Real),Dict{Int64,Float64}}, Tuple{AbstractArray{#s201,2} where #s201<:Real,AbstractArray{#s200,1} where #s200<:(AbstractArray{#s199,1} where #s199<:Real),Dict{Int64,Float64},Tuple{Vararg{Int64,N}}}} where N",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.intralinks",
    "category": "method",
    "text": "intralinks( cluster::AbstractMatrix{<:Real},\n            vectors::AbstractVector{<:AbstractVector{<:Real}},\n            neighbors::Dict{Int,Float},\n            maxtranslations::NTuple{N,Int}=ntuple(i->length(neighbors),length(vectors))\n            ) where N -> Vector{Link}\n\nUse kdtree to get the intracluster nearest neighbors.\n\nAs is similar to minimumlengths, when vectors is nonempty, the cluster assumes periodic boundaries. neighbors provides the map between the bond length and the order of nearest neighbors. Note only those with the lengths present in neighbors will be included in the result. maxtranslations determines the maximum number of translations along those directions specified by vectors when the tiled supercluster is construted (See minimumlengths for the explanation of the method for periodic lattices).\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.isintracell-Tuple{Bond}",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.isintracell",
    "category": "method",
    "text": "isintracell(bond::Bond) -> Bool\n\nJudge whether a bond is intra the unit cell of a lattice.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.isintratriangle-Tuple{AbstractArray{#s15,1} where #s15<:Real,AbstractArray{#s14,1} where #s14<:Real,AbstractArray{#s13,1} where #s13<:Real,AbstractArray{#s12,1} where #s12<:Real}",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.isintratriangle",
    "category": "method",
    "text": "isintratriangle(p::AbstractVector{<:Real},\n                p1::AbstractVector{<:Real},\n                p2::AbstractVector{<:Real},\n                p3::AbstractVector{<:Real};\n                vertexes::NTuple{3,Bool}=(true,true,true),\n                edges::NTuple{3,Bool}=(true,true,true),\n                atol::Real=atol,\n                rtol::Real=rtol\n                ) -> Bool\n\nJudge whether a point belongs to the interior of a triangle whose vertexes are p1, \'p2\' and p3 with the give tolerance. vertexes and edges define whether the interior should contain the vertexes or edges, respectively.\n\nnotes: Notes\nThe vertexes are in the order (p1,p2,p3) and the edges are in the order (p1p2,p2p3,p3p1).\nThe edges do not contain the vertexes.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.isonline-Tuple{AbstractArray{#s15,1} where #s15<:Real,AbstractArray{#s14,1} where #s14<:Real,AbstractArray{#s13,1} where #s13<:Real}",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.isonline",
    "category": "method",
    "text": "isonline(p::AbstractVector{<:Real},p1::AbstractVector{<:Real},p2::AbstractVector{<:Real};ends::Tuple{Bool,Bool}=(true,true),atol::Real=atol,rtol::Real=rtol) -> Bool\n\nJudge whether a point is on a line segment whose end points are p1 and p2 with the given tolerance. ends defines whether the line segment should contain its ends.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.isparallel-Tuple{AbstractArray{#s68,1} where #s68<:Real,AbstractArray{#s67,1} where #s67<:Real}",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.isparallel",
    "category": "method",
    "text": "isparallel(v1::AbstractVector{<:Real},v2::AbstractVector{<:Real};atol::Real=atol,rtol::Real=rtol) -> Bool\n\nJudge whether two vectors are parallel to each other with the given tolerance, 0 for not parallel, 1 for parallel and -1 for antiparallel.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.issubordinate-Tuple{AbstractArray{#s203,1} where #s203<:Real,AbstractArray{#s202,1} where #s202<:(AbstractArray{#s201,1} where #s201<:Real)}",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.issubordinate",
    "category": "method",
    "text": "issubordinate(rcoord::AbstractVector{<:Real},vectors::AbstractVector{<:AbstractVector{<:Real}};atol::Real=atol,rtol::Real=rtol) -> Bool\n\nJudge whether a coordinate belongs to a lattice defined by vectors with the given tolerance.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.minimumlengths",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.minimumlengths",
    "category": "function",
    "text": "minimumlengths(cluster::AbstractMatrix{<:Real},vectors::AbstractVector{<:AbstractVector{<:Real}},nneighbor::Int=1;coordination::Int=8) -> Vector{Float}\n\nUse kdtree to search the lowest several minimum bond lengths within a lattice translated by a cluster.\n\nWhen the translation vectors are not empty, the lattice will be considered periodic in the corresponding directions. Otherwise the lattice will be open in all directions. To search for the bonds accorss the periodic boundaries, the cluster will be pretranslated to become a supercluster, which has open boundaries but is large enough to contain all the nearest neighbors within the required order. The coordination parameter sets the average number of each order of nearest neighbors. If it is to small, larger bond lengths may not be searched, and the result will contain Inf. This is a sign that you may need a larger coordination. Another situation that Inf appears in the result occurs when the minimum lengths are searched in open lattices. Indeed, the cluster may be too small so that the required order just goes beyond it. In this case the warning message can be safely ignored.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.nneighbor-Tuple{AbstractLattice}",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.nneighbor",
    "category": "method",
    "text": "nneighbor(lattice::AbstractLattice) -> Int\n\nGet the highest order of nearest neighbors.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.pidtype-Tuple{AbstractBond}",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.pidtype",
    "category": "method",
    "text": "pidtype(b::AbstractBond)\npidtype(::Type{<:AbstractBond{R,P,N}}) where {R,P,N}\n\nGet the pid type of a concrete bond.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.polar-Tuple{AbstractArray{#s199,1} where #s199<:Real}",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.polar",
    "category": "method",
    "text": "polar(v::AbstractVector{<:Real}) -> Float\n\nGet the polar angle in radians of a vector.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.polard-Tuple{AbstractArray{#s199,1} where #s199<:Real}",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.polard",
    "category": "method",
    "text": "polard(v::AbstractVector{<:Real}) -> Float\n\nGet the polar angle in degrees of a vector.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.rcoord-Tuple{Bond}",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.rcoord",
    "category": "method",
    "text": "rcoord(bond::Bond) -> SVector\n\nGet the rcoord of the bond.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.reciprocals-Tuple{AbstractArray{#s213,1} where #s213<:(AbstractArray{#s212,1} where #s212<:Real)}",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.reciprocals",
    "category": "method",
    "text": "reciprocals(vectors::AbstractVector{AbstractVector{<:Real}}) -> Vector{Vector{Float}}\n\nGet the reciprocals dual to the input vectors.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.rotate-Tuple{AbstractArray{#s202,2} where #s202<:Real,Real}",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.rotate",
    "category": "method",
    "text": "rotate(cluster::AbstractMatrix{<:Real},angle::Real;axis::Tuple{Union{AbstractVector{<:Real},Nothing},Tuple{<:Real,<:Real}}=(nothing,(0,0))) -> Matrix{Float}\n\nGet a rotated cluster of the original one by a certain angle around an axis.\n\nThe axis is determined by a point it gets through (nothing can be used to denote the origin), and its polar as well as azimuth angles in radians. The default axis is the z axis.\n\nnotes: Notes\nThe result is given by the Rodrigues\' rotation formula.\nOnly 2 and 3 dimensional vectors can be rotated.\nWhen the input vectors are 2 dimensional, both the polar and azimuth of the axis must be 0.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.tile-Union{Tuple{M}, Tuple{N}, Tuple{AbstractArray{#s204,2} where #s204<:Real,AbstractArray{#s203,1} where #s203<:(AbstractArray{#s202,1} where #s202<:Real)}, Tuple{AbstractArray{#s201,2} where #s201<:Real,AbstractArray{#s200,1} where #s200<:(AbstractArray{#s199,1} where #s199<:Real),Tuple{Vararg{Tuple{Vararg{#s198,N}} where #s198<:Real,M}}}} where M where N",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.tile",
    "category": "method",
    "text": "tile(cluster::AbstractMatrix{<:Real},vectors::AbstractVector{<:AbstractVector{<:Real}},translations::NTuple{M,NTuple{N,<:Real}}=()) where {N,M} -> Matrix{Float}\n\nTile a supercluster by translations of the input cluster.\n\nBasically, the final supercluster is composed of several parts, each of which is a translation of the original cluster, with the translation vectors specified by vectors and each set of the translation indices contained in translations. When translation vectors are empty, a copy of the original cluster will be returned.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.translate-Tuple{AbstractArray{#s213,2} where #s213<:Real,AbstractArray{#s212,1} where #s212<:Real}",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.translate",
    "category": "method",
    "text": "translate(cluster::AbstractMatrix{<:Real},vector::AbstractVector{<:Real}) -> Matrix{vector|>eltype}\n\nGet the translated cluster of the original one by a vector.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Essentials.Spatial.volume-Tuple{AbstractArray{#s188,1} where #s188<:Real,AbstractArray{#s187,1} where #s187<:Real,AbstractArray{#s186,1} where #s186<:Real}",
    "page": "Spatial",
    "title": "Hamiltonian.Essentials.Spatial.volume",
    "category": "method",
    "text": "volume(v1::AbstractVector{<:Real},v2::AbstractVector{<:Real},v3::AbstractVector{<:Real}) -> Real\n\nGet the volume spanned by three vectors.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Utilities.dimension-Tuple{AbstractBond}",
    "page": "Spatial",
    "title": "Hamiltonian.Utilities.dimension",
    "category": "method",
    "text": "dimension(b::AbstractBond) -> Int\ndimension(::Type{<:AbstractBond{R,P,N}}) where {R,P,N} -> Int\n\nGet the space dimension of a concrete bond.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Utilities.dimension-Tuple{AbstractLattice}",
    "page": "Spatial",
    "title": "Hamiltonian.Utilities.dimension",
    "category": "method",
    "text": "dimension(lattice::AbstractLattice) -> Int\ndimension(::Type{<:AbstractLattice{P,N}}) where {P,N}-> Int\n\nGet the space dimension of the lattice.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Utilities.dimension-Tuple{Link}",
    "page": "Spatial",
    "title": "Hamiltonian.Utilities.dimension",
    "category": "method",
    "text": "dimension(link::Link) -> Int\ndimension(::Type{<:Link{N}}) where N -> Int\n\nGet the space dimension of a link.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Hamiltonian.Utilities.rank-Tuple{AbstractBond}",
    "page": "Spatial",
    "title": "Hamiltonian.Utilities.rank",
    "category": "method",
    "text": "rank(b::AbstractBond) -> Int\nrank(::Type{<:AbstractBond{R,P,N}}) where {R,P,N} -> Int\n\nGet the rank of a concrete bond.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Base.:==-Tuple{AbstractLattice,AbstractLattice}",
    "page": "Spatial",
    "title": "Base.:==",
    "category": "method",
    "text": "==(lattice1::AbstractLattice,lattice2::AbstractLattice) -> Bool\nisequal(lattice1::AbstractLattice,lattice2::AbstractLattice) -> Bool\n\nOverloaded equivalent operator.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Base.:==-Tuple{Link,Link}",
    "page": "Spatial",
    "title": "Base.:==",
    "category": "method",
    "text": "==(l1::Link,l2::Link) -> Bool\nisequal(l1::Link,l2::Link) -> Bool\n\nOverloaded equivalent operator.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Base.:==-Union{Tuple{R}, Tuple{AbstractBond{R,P,N} where N where P<:PID,AbstractBond{R,P,N} where N where P<:PID}} where R",
    "page": "Spatial",
    "title": "Base.:==",
    "category": "method",
    "text": "==(b1::AbstractBond{R},b1::AbstractBond{R}) where R -> Bool\nisequal(b1::AbstractBond{R},b1::AbstractBond{R}) where R -> Bool\n\nOverloaded equivalent operator.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Base.getindex-Tuple{AbstractLattice,RCoordIndex{Int64}}",
    "page": "Spatial",
    "title": "Base.getindex",
    "category": "method",
    "text": "getindex(lattice::AbstractLattice,i::RCoordIndex) -> SVector\ngetindex(lattice::AbstractLattice,i::ICoordIndex) -> SVector\ngetindex(lattice::AbstractLattice,i::PointIndex) -> Point\n\nGet a rcoord, an icoord or a point of a lattice according to the type of the input index.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Base.insert!-Union{Tuple{S}, Tuple{Cylinder,Vararg{S,N} where N}} where S",
    "page": "Spatial",
    "title": "Base.insert!",
    "category": "method",
    "text": "insert!(cylinder::Cylinder,ps::S...;cut::Int=length(cylinder)÷2+1,scopes::Union{<:AbstractVector{S},Nothing}=nothing,coordination::Int=9) where S -> Cylinder\n\nInsert a couple of blocks into a cylinder.\n\nThe position of the cut of the cylinder is specified by the keyword argument cut, which is the center of the cylinder by default. All pids corresponding to a same newly inserted block share the same scope, which is specified by the parameter ps. Optionally, the scopes of the old pids in the cylinder can be replaced if the parameter scopes is assigned other than nothing. Note the length of ps is equal to the number of newly inserted blocks, while that of scopes should be equal to the old length of the cylinder.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Base.keytype-Tuple{AbstractLattice}",
    "page": "Spatial",
    "title": "Base.keytype",
    "category": "method",
    "text": "keytype(lattice::AbstractLattice)\nkeytype(::Type{<:AbstractLattice{P,N}}) where {P,N}\n\nGet the pid type of the lattice.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Base.length-Tuple{AbstractLattice}",
    "page": "Spatial",
    "title": "Base.length",
    "category": "method",
    "text": "length(lattice::AbstractLattice) -> Int\n\nGet the number of points contained in a lattice.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Base.reverse-Tuple{Bond}",
    "page": "Spatial",
    "title": "Base.reverse",
    "category": "method",
    "text": "reverse(bond::Bond) -> Bond\n\nGet the reversed bond.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Base.show-Tuple{IO,AbstractLattice}",
    "page": "Spatial",
    "title": "Base.show",
    "category": "method",
    "text": "show(io::IO,lattice::AbstractLattice)\n\nShow a lattice.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Base.show-Tuple{IO,Bond}",
    "page": "Spatial",
    "title": "Base.show",
    "category": "method",
    "text": "show(io::IO,bond::Bond)\n\nShow a bond.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Base.show-Tuple{IO,Link}",
    "page": "Spatial",
    "title": "Base.show",
    "category": "method",
    "text": "show(io::IO,link::Link)\n\nShow a link.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Base.show-Tuple{IO,Point}",
    "page": "Spatial",
    "title": "Base.show",
    "category": "method",
    "text": "show(io::IO,p::Point)\n\nShow a labeled point.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Base.valtype-Tuple{AbstractLattice}",
    "page": "Spatial",
    "title": "Base.valtype",
    "category": "method",
    "text": "valtype(lattice::AbstractLattice)\nvaltype(::Type{<:AbstractLattice{P,N}}) where {P,N}\n\nGet the point type of the lattice.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatial.html#Manul-1",
    "page": "Spatial",
    "title": "Manul",
    "category": "section",
    "text": "Modules=[Spatial]\nOrder=  [:module,:constant,:type,:macro,:function]"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#",
    "page": "DegreeOfFreedom",
    "title": "DegreeOfFreedom",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Essentials.DegreeOfFreedompush!(LOAD_PATH,\"../../../../src/\")\nusing Hamiltonian"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#DegreeOfFreedom-1",
    "page": "DegreeOfFreedom",
    "title": "DegreeOfFreedom",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Table-1",
    "page": "DegreeOfFreedom",
    "title": "Table",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Internal-and-Index-1",
    "page": "DegreeOfFreedom",
    "title": "Internal and Index",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#IDFConfig-1",
    "page": "DegreeOfFreedom",
    "title": "IDFConfig",
    "category": "section",
    "text": "There are two purposes with IDFConfigprovide a complete set of internal degrees of freedom on a lattice\noffer the ordering of the internal degrees of freedom on a lattice"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#IndexPack-and-IndexPacks-1",
    "page": "DegreeOfFreedom",
    "title": "IndexPack and IndexPacks",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Hamiltonian.Essentials.DegreeOfFreedom.Coupling",
    "page": "DegreeOfFreedom",
    "title": "Hamiltonian.Essentials.DegreeOfFreedom.Coupling",
    "category": "type",
    "text": "Coupling{V,I,N} <: Element{V,I,N}\n\nThe coupling intra/inter interanl degrees of freedom at different lattice points.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Hamiltonian.Essentials.DegreeOfFreedom.Couplings",
    "page": "DegreeOfFreedom",
    "title": "Hamiltonian.Essentials.DegreeOfFreedom.Couplings",
    "category": "type",
    "text": "Couplings{I<:ID,C<:Coupling} <: AbstractDict{I,C}\n\nA pack of couplings intra/inter interanl degrees of freedom at different lattice points.\n\nAlias for Elements{I,C}.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Hamiltonian.Essentials.DegreeOfFreedom.DirectIndexToTuple",
    "page": "DegreeOfFreedom",
    "title": "Hamiltonian.Essentials.DegreeOfFreedom.DirectIndexToTuple",
    "category": "type",
    "text": "DirectIndexToTuple\n\nDirect index to tuple.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Hamiltonian.Essentials.DegreeOfFreedom.DirectIndexToTuple-Tuple{Index}",
    "page": "DegreeOfFreedom",
    "title": "Hamiltonian.Essentials.DegreeOfFreedom.DirectIndexToTuple",
    "category": "method",
    "text": "(indextotuple::DirectIndexToTuple)(index::Index) -> Tuple\n\nConvert an index to tuple directly.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Hamiltonian.Essentials.DegreeOfFreedom.FilteredAttributes",
    "page": "DegreeOfFreedom",
    "title": "Hamiltonian.Essentials.DegreeOfFreedom.FilteredAttributes",
    "category": "type",
    "text": "FilteredAttributes(::Type{I}) where I<:Index\n\nA method that converts an arbitary index to a tuple, by iterating over the selected attributes in a specific order.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Hamiltonian.Essentials.DegreeOfFreedom.FilteredAttributes-Tuple{Index}",
    "page": "DegreeOfFreedom",
    "title": "Hamiltonian.Essentials.DegreeOfFreedom.FilteredAttributes",
    "category": "method",
    "text": "(indextotuple::FilteredAttributes)(index::Index) -> Tuple\n\nConvert an index to tuple by the \"filtered attributes\" method.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Hamiltonian.Essentials.DegreeOfFreedom.IDFConfig",
    "page": "DegreeOfFreedom",
    "title": "Hamiltonian.Essentials.DegreeOfFreedom.IDFConfig",
    "category": "type",
    "text": "IDFConfig(map::Function,::Type{I},pids::AbstractVector{<:PID}=[]) where I<:Internal\n\nConfiguration of the internal degrees of freedom at a lattice.\n\nmap maps a PID to an Internal.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Hamiltonian.Essentials.DegreeOfFreedom.IID",
    "page": "DegreeOfFreedom",
    "title": "Hamiltonian.Essentials.DegreeOfFreedom.IID",
    "category": "type",
    "text": "IID\n\nThe id of an internal degree of freedom.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Hamiltonian.Essentials.DegreeOfFreedom.Index",
    "page": "DegreeOfFreedom",
    "title": "Hamiltonian.Essentials.DegreeOfFreedom.Index",
    "category": "type",
    "text": "Index{P,I}\n\nThe complete index of a degree of freedom, which consist of the spatial part and the internal part.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Hamiltonian.Essentials.DegreeOfFreedom.IndexToTuple",
    "page": "DegreeOfFreedom",
    "title": "Hamiltonian.Essentials.DegreeOfFreedom.IndexToTuple",
    "category": "type",
    "text": "IndexToTuple\n\nThe rules for converting an index to a tuple.\n\nAs a function, every instance should accept only one positional argument, i.e. the index to be converted to a tuple.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Hamiltonian.Essentials.DegreeOfFreedom.Internal",
    "page": "DegreeOfFreedom",
    "title": "Hamiltonian.Essentials.DegreeOfFreedom.Internal",
    "category": "type",
    "text": "Internal\n\nThe whole internal degrees of freedom at a single point.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Hamiltonian.Essentials.DegreeOfFreedom.Table",
    "page": "DegreeOfFreedom",
    "title": "Hamiltonian.Essentials.DegreeOfFreedom.Table",
    "category": "type",
    "text": "Table{I<:Index} <: AbstractDict{I,Int}\n\nIndex-sequence table. Alias for Dict{I<:Index,Int}.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Hamiltonian.Essentials.DegreeOfFreedom.Table-Tuple{AbstractArray{#s211,1} where #s211<:Index}",
    "page": "DegreeOfFreedom",
    "title": "Hamiltonian.Essentials.DegreeOfFreedom.Table",
    "category": "method",
    "text": "Table(indices::AbstractVector{<:Index};by::IndexToTuple=directindextotuple) -> Table\n\nConvert an sequence of indices to the corresponding index-sequence table.\n\nThe input indices will be converted to tuples by the by function with the duplicates removed. The resulting unique tuples are sorted, which determines the sequence of the input indices. Note that two indices have the same sequence if their converted tupels are equal to each other.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Hamiltonian.Essentials.DegreeOfFreedom.Table-Tuple{IDFConfig}",
    "page": "DegreeOfFreedom",
    "title": "Hamiltonian.Essentials.DegreeOfFreedom.Table",
    "category": "method",
    "text": "Table(config::IDFConfig;by::IndexToTuple=directindextotuple) -> Table\n\nGet the index-sequence table of the whole internal Hilbert spaces at a lattice.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Hamiltonian.Essentials.DegreeOfFreedom.directindextotuple",
    "page": "DegreeOfFreedom",
    "title": "Hamiltonian.Essentials.DegreeOfFreedom.directindextotuple",
    "category": "function",
    "text": "directindextotuple\n\nIndicate that the conversion from an index to a tuple is direct.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Hamiltonian.Essentials.DegreeOfFreedom.iid-Tuple{Index}",
    "page": "DegreeOfFreedom",
    "title": "Hamiltonian.Essentials.DegreeOfFreedom.iid",
    "category": "method",
    "text": "iid(index::Index) -> IID\n\nGet the internal part of an index.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Hamiltonian.Essentials.DegreeOfFreedom.iidtype-Tuple{Index}",
    "page": "DegreeOfFreedom",
    "title": "Hamiltonian.Essentials.DegreeOfFreedom.iidtype",
    "category": "method",
    "text": "iidtype(index::Index)\niidtype(::Type{<:Index{P,I}}) where {P,I}\n\nGet the type of the internal part of an index.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Hamiltonian.Essentials.DegreeOfFreedom.pid-Tuple{Index}",
    "page": "DegreeOfFreedom",
    "title": "Hamiltonian.Essentials.DegreeOfFreedom.pid",
    "category": "method",
    "text": "pid(index::Index) -> PID\n\nGet the spatial part of an index.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Hamiltonian.Essentials.Spatial.pidtype-Tuple{Index}",
    "page": "DegreeOfFreedom",
    "title": "Hamiltonian.Essentials.Spatial.pidtype",
    "category": "method",
    "text": "pidtype(index::Index)\npidtype(::Type{<:Index{P,I}}) where {P,I}\n\nGet the type of the spatial part of an index.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Core.Type-Tuple{PID,IID}",
    "page": "DegreeOfFreedom",
    "title": "Core.Type",
    "category": "method",
    "text": "(INDEX::Type{<:Index})(pid::PID,iid::IID) -> INDEX\n\nGet the corresponding index from a pid and an iid.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Base.:==-Union{Tuple{I}, Tuple{I,I}} where I<:Internal",
    "page": "DegreeOfFreedom",
    "title": "Base.:==",
    "category": "method",
    "text": "==(i1::I,i2::I) where I<:Internal -> Bool\nisequal(i1::I,i2::I) where I<:Internal -> Bool\n\nCompare two internals and judge whether they are equal to each other.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Base.adjoint-Tuple{Index}",
    "page": "DegreeOfFreedom",
    "title": "Base.adjoint",
    "category": "method",
    "text": "adjoint(index::Index) -> typeof(index)\n\nGet the adjoint of an index.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Base.eltype-Tuple{Internal}",
    "page": "DegreeOfFreedom",
    "title": "Base.eltype",
    "category": "method",
    "text": "eltype(internal::Internal)\neltype(::Type{<:Internal{I}}) where I\n\nGet the type of the IIDs that an Internal contains.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Base.filter-Tuple{Function,FilteredAttributes}",
    "page": "DegreeOfFreedom",
    "title": "Base.filter",
    "category": "method",
    "text": "filter(f::Function,indextotuple::FilteredAttributes)\n\nFilter the attributes of a \"filtered attributes\" method.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Base.length-Tuple{FilteredAttributes}",
    "page": "DegreeOfFreedom",
    "title": "Base.length",
    "category": "method",
    "text": "length(indextotuple::FilteredAttributes) -> Int\nlength(::Type{<:FilteredAttributes{N}}) where N -> Int\n\nGet the length of the filtered attributes.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Base.replace!-Tuple{IDFConfig,Vararg{PID,N} where N}",
    "page": "DegreeOfFreedom",
    "title": "Base.replace!",
    "category": "method",
    "text": "replace!(config::IDFConfig,pids::PID...) -> IDFConfig\n\nReset the idfconfig with new pids.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Base.reverse-Tuple{Dict{I,Int64} where I<:Index}",
    "page": "DegreeOfFreedom",
    "title": "Base.reverse",
    "category": "method",
    "text": "reverse(table::Table) -> Dict{Int,Set{<:Index}}\n\nConvert an index-sequence table to a sequence-indices table.\n\nSince different indices may correspond to the same sequence, the reverse is a one-to-many map.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Base.show-Tuple{IO,Internal}",
    "page": "DegreeOfFreedom",
    "title": "Base.show",
    "category": "method",
    "text": "show(io::IO,i::Internal)\n\nShow an internal.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Base.union-Tuple{Vararg{Dict{I,Int64} where I<:Index,N} where N}",
    "page": "DegreeOfFreedom",
    "title": "Base.union",
    "category": "method",
    "text": "union(tables::Table...;by::IndexToTuple=directindextotuple) -> Table\n\nUnite several index-sequence tables.\n\nSee Table for more details.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Base.union-Union{Tuple{I}, Tuple{P}, Tuple{Type{P},Type{I}}} where I<:IID where P<:PID",
    "page": "DegreeOfFreedom",
    "title": "Base.union",
    "category": "method",
    "text": "union(::Type{P},::Type{I}) where {P<:PID,I<:IID}\n\nCombine a concrete PID type and a concrete IID type to a concrete Index type.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreeOfFreedom.html#Manul-1",
    "page": "DegreeOfFreedom",
    "title": "Manul",
    "category": "section",
    "text": "Modules=[DegreeOfFreedom]\nOrder=  [:module,:constant,:type,:macro,:function]"
},

{
    "location": "man/Essentials/FockPackage.html#",
    "page": "Fock Package",
    "title": "Fock Package",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Essentials.FockPackagepush!(LOAD_PATH,\"../../../../src/\")\nusing Hamiltonian"
},

{
    "location": "man/Essentials/FockPackage.html#Fock-Package-1",
    "page": "Fock Package",
    "title": "Fock Package",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.ANNIHILATION",
    "page": "Fock Package",
    "title": "Hamiltonian.Essentials.FockPackage.ANNIHILATION",
    "category": "constant",
    "text": "ANNIHILATION\n\nIndicate that the nambu index is ANNIHILATION.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.CREATION",
    "page": "Fock Package",
    "title": "Hamiltonian.Essentials.FockPackage.CREATION",
    "category": "constant",
    "text": "CREATION\n\nIndicate that the nambu index is CREATION.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.FCID",
    "page": "Fock Package",
    "title": "Hamiltonian.Essentials.FockPackage.FCID",
    "category": "type",
    "text": "FCID{A,O,S,N} <: SimpleID\n\nThe id of a Fock coupling.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.FCID-Tuple{}",
    "page": "Fock Package",
    "title": "Hamiltonian.Essentials.FockPackage.FCID",
    "category": "method",
    "text": "FCID(;atom=nothing,orbital=nothing,spin=nothing,nambu=nothing)\n\nConstruct an id of a Fock coupling.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.FID",
    "page": "Fock Package",
    "title": "Hamiltonian.Essentials.FockPackage.FID",
    "category": "type",
    "text": "FID <: IID\n\nThe Fock id.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.FID-Tuple{}",
    "page": "Fock Package",
    "title": "Hamiltonian.Essentials.FockPackage.FID",
    "category": "method",
    "text": "FID(;orbital::Int=1,spin::Int=1,nambu::Int=ANNIHILATION)\n\nCreate a Fock id.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.FIndex",
    "page": "Fock Package",
    "title": "Hamiltonian.Essentials.FockPackage.FIndex",
    "category": "type",
    "text": "FIndex{S} <: Index{PID{S},FID}\n\nThe Fock index.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.Fock",
    "page": "Fock Package",
    "title": "Hamiltonian.Essentials.FockPackage.Fock",
    "category": "type",
    "text": "Fock <: Internal{FID}\n\nThe Fock interanl degrees of freedom.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.Fock-Tuple{}",
    "page": "Fock Package",
    "title": "Hamiltonian.Essentials.FockPackage.Fock",
    "category": "method",
    "text": "Fock(;atom::Int=1,norbital::Int=1,nspin::Int=2,nnambu::Int=2)\n\nConstruct a Fock degrees of freedom.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.FockCoupling",
    "page": "Fock Package",
    "title": "Hamiltonian.Essentials.FockPackage.FockCoupling",
    "category": "type",
    "text": "FockCoupling{N,V,I<:ID{N,<:FCID}} <: Coupling{V,I,N}\n\nA Fock coupling.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.FockCoupling-Union{Tuple{Number}, Tuple{N}} where N",
    "page": "Fock Package",
    "title": "Hamiltonian.Essentials.FockPackage.FockCoupling",
    "category": "method",
    "text": "FockCoupling{N}(    value::Number;\n                    atoms::Union{NTuple{N,Int},Nothing}=nothing,\n                    orbitals::Union{NTuple{N,Int},Nothing}=nothing,\n                    spins::Union{NTuple{N,Int},Nothing}=nothing,\n                    nambus::Union{NTuple{N,Int},Nothing}=nothing) where N\n\nConstruct a Fock coupling.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.nambufockindextotuple",
    "page": "Fock Package",
    "title": "Hamiltonian.Essentials.FockPackage.nambufockindextotuple",
    "category": "function",
    "text": "nambufockindextotuple\n\nIndicate that the filtered attributes are (:scope,:nambu,:site,:orbital,:spin) when converting a Fock index to tuple.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.usualfockindextotuple",
    "page": "Fock Package",
    "title": "Hamiltonian.Essentials.FockPackage.usualfockindextotuple",
    "category": "function",
    "text": "usualfockindextotuple\n\nIndicate that the filtered attributes are (:scope,:site,:orbital,:spin) when converting a Fock index to tuple.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.σʸ-Tuple{String}",
    "page": "Fock Package",
    "title": "Hamiltonian.Essentials.FockPackage.σʸ",
    "category": "method",
    "text": "σʸ(mode::String)\n\nThe Pauli matrix σʸ, which can act on the space of spins(\"sp\"), orbitals(\"ob\"), sublattices(\"sl\") or particle-holes(\"ph\").\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.σˣ-Tuple{String}",
    "page": "Fock Package",
    "title": "Hamiltonian.Essentials.FockPackage.σˣ",
    "category": "method",
    "text": "σˣ(mode::String)\n\nThe Pauli matrix σˣ, which can act on the space of spins(\"sp\"), orbitals(\"ob\"), sublattices(\"sl\") or particle-holes(\"ph\").\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.σᶻ-Tuple{String}",
    "page": "Fock Package",
    "title": "Hamiltonian.Essentials.FockPackage.σᶻ",
    "category": "method",
    "text": "σᶻ(mode::String)\n\nThe Pauli matrix σᶻ, which can act on the space of spins(\"sp\"), orbitals(\"ob\"), sublattices(\"sl\") or particle-holes(\"ph\").\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.σ⁰-Tuple{String}",
    "page": "Fock Package",
    "title": "Hamiltonian.Essentials.FockPackage.σ⁰",
    "category": "method",
    "text": "σ⁰(mode::String)\n\nThe Pauli matrix σ⁰, which can act on the space of spins(\"sp\"), orbitals(\"ob\"), sublattices(\"sl\") or particle-holes(\"ph\").\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.σ⁺-Tuple{String}",
    "page": "Fock Package",
    "title": "Hamiltonian.Essentials.FockPackage.σ⁺",
    "category": "method",
    "text": "σ⁺(mode::String)\n\nThe Pauli matrix σ⁺, which can act on the space of spins(\"sp\"), orbitals(\"ob\"), sublattices(\"sl\") or particle-holes(\"ph\").\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.σ⁻-Tuple{String}",
    "page": "Fock Package",
    "title": "Hamiltonian.Essentials.FockPackage.σ⁻",
    "category": "method",
    "text": "σ⁻(mode::String)\n\nThe Pauli matrix σ⁻, which can act on the space of spins(\"sp\"), orbitals(\"ob\"), sublattices(\"sl\") or particle-holes(\"ph\").\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Utilities.expand-Tuple{FockCoupling,AbstractBond,Vararg{Fock,N} where N}",
    "page": "Fock Package",
    "title": "Hamiltonian.Utilities.expand",
    "category": "method",
    "text": "expand(fc::FockCoupling,bond::AbstractBond,focks::Fock...) -> FCExpand\n\nExpand a Fock coupling on a bond with the given Fock degrees of freedom.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Base.:*-Tuple{FockCoupling{2,V,I} where I<:(Hamiltonian.Utilities.AlgebraOverField.ID{2,#s232} where #s232<:FCID) where V,FockCoupling{2,V,I} where I<:(Hamiltonian.Utilities.AlgebraOverField.ID{2,#s232} where #s232<:FCID) where V}",
    "page": "Fock Package",
    "title": "Base.:*",
    "category": "method",
    "text": "*(fc1::FockCoupling{2},fc2::FockCoupling{2}) -> FockCoupling{2}\n\nGet the multiplication between two rank-2 Fock couplings.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Base.adjoint-Tuple{FID}",
    "page": "Fock Package",
    "title": "Base.adjoint",
    "category": "method",
    "text": "adjoint(fid::FID) -> FID\n\nGet the adjoint of a Fock id.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Base.fieldnames-Tuple{Type{#s215} where #s215<:FCID}",
    "page": "Fock Package",
    "title": "Base.fieldnames",
    "category": "method",
    "text": "fieldnames(::Type{<:FCID}) -> NTuple{4,Symbol}\n\nGet the fieldnames of FCID.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Base.fieldnames-Tuple{Type{#s215} where #s215<:FIndex}",
    "page": "Fock Package",
    "title": "Base.fieldnames",
    "category": "method",
    "text": "fieldnames(::Type{<:FIndex}) -> NTuple{5,Symbol}\n\nGet the fieldnames of a Fock index.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Base.fieldnames-Tuple{Type{FID}}",
    "page": "Fock Package",
    "title": "Base.fieldnames",
    "category": "method",
    "text": "fieldnames(::Type{FID}) -> NTuple{3,Symbol}\n\nGet the fieldnames of FID.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Base.iterate",
    "page": "Fock Package",
    "title": "Base.iterate",
    "category": "function",
    "text": "iterate(fock::Fock,state=1)\n\nIterate over a Fock degrees of freedom.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Base.length-Tuple{Fock}",
    "page": "Fock Package",
    "title": "Base.length",
    "category": "method",
    "text": "length(fock::Fock) -> Int\n\nGet the dimension of a Fock degrees of freedom.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Base.repr-Tuple{FockCoupling}",
    "page": "Fock Package",
    "title": "Base.repr",
    "category": "method",
    "text": "repr(fc::FockCoupling) -> String\n\nGet the repr representation of a Fock coupling.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Base.show-Tuple{IO,FockCoupling}",
    "page": "Fock Package",
    "title": "Base.show",
    "category": "method",
    "text": "show(io::IO,fc::FockCoupling)\n\nShow a Fock coupling.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Base.union-Union{Tuple{P}, Tuple{Type{P},Type{FID}}} where P<:PID",
    "page": "Fock Package",
    "title": "Base.union",
    "category": "method",
    "text": "union(::Type{P},::Type{FID})\n\nGet the union type of PID and FID.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Fock-Degree-of-Freedom-1",
    "page": "Fock Package",
    "title": "Fock Degree of Freedom",
    "category": "section",
    "text": "Modules=[FockPackage]\nOrder=  [:module,:constant,:type,:macro,:function]"
},

]}
