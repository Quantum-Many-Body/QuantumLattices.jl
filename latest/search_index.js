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
    "location": "man/Prerequisites/Introduction.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Prerequisites"
},

{
    "location": "man/Prerequisites/Introduction.html#Introduction-1",
    "page": "Introduction",
    "title": "Introduction",
    "category": "section",
    "text": "This module contains the prerequisites of the Hamiltonian package.The constants, types, macros, functions or submodules defined in this module will NOT be exported by the package. Instead, they serve as the prerequisites and fundamentals. The range of the contents are quite wide, but basically, they fall into 3 categories:Global constants and miscellaneous tiny useful functions;\nGeneric functions that are extended by other parts of the package;\nBasic data structures as supplements to the Julia.Base and other common packages.The first category is contained in the main body of this module, while the others come in separate submodules."
},

{
    "location": "man/Prerequisites/Introduction.html#Hamiltonian.Prerequisites.atol",
    "page": "Introduction",
    "title": "Hamiltonian.Prerequisites.atol",
    "category": "constant",
    "text": "Absolute tolerance for float numbers.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Introduction.html#Hamiltonian.Prerequisites.rtol",
    "page": "Introduction",
    "title": "Hamiltonian.Prerequisites.rtol",
    "category": "constant",
    "text": "Relative tolerance for float numbers.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Introduction.html#Hamiltonian.Prerequisites.Float",
    "page": "Introduction",
    "title": "Hamiltonian.Prerequisites.Float",
    "category": "type",
    "text": "Default float type.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Introduction.html#Hamiltonian.Prerequisites.decimaltostr",
    "page": "Introduction",
    "title": "Hamiltonian.Prerequisites.decimaltostr",
    "category": "function",
    "text": "decimaltostr(number::Integer,n::Int=5)\ndecimaltostr(number::Rational,n::Int=5)\ndecimaltostr(number::AbstractFloat,n::Int=5)\ndecimaltostr(number::Complex,n::Int=5)\n\nConvert a number to a string with at most n decimal places.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Introduction.html#Hamiltonian.Prerequisites.ordinal",
    "page": "Introduction",
    "title": "Hamiltonian.Prerequisites.ordinal",
    "category": "function",
    "text": "ordinal(number::Interger)\n\nConvert a positive number to its corresponding ordinal.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Introduction.html#Hamiltonian.Prerequisites.delta",
    "page": "Introduction",
    "title": "Hamiltonian.Prerequisites.delta",
    "category": "function",
    "text": "delta(i,j) -> Int\n\nKronecker delta function.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Introduction.html#Constants-and-functions-1",
    "page": "Introduction",
    "title": "Constants and functions",
    "category": "section",
    "text": "All the following constants and functions in this section are defined in the main body and are exported by this module.atol\nrtol\nFloat\ndecimaltostr\nordinal\ndelta"
},

{
    "location": "man/Prerequisites/Introduction.html#Generic-interfaces-1",
    "page": "Introduction",
    "title": "Generic interfaces",
    "category": "section",
    "text": "For details, please refer to the page Interfaces:Pages=[\n    \"Interfaces.md\",\n    ]\nDepth=2"
},

{
    "location": "man/Prerequisites/Introduction.html#Basic-structures-1",
    "page": "Introduction",
    "title": "Basic structures",
    "category": "section",
    "text": "Here lists the table of contents of the basic data structures that are supplements to the Julia.Base and other common packages:Pages=[\n    \"TypeTraits.md\",\n    \"Factories.md\",\n    \"CompositeStructures.md\",\n    \"Trees.md\",\n    \"NamedVectors.md\",\n    ]\nDepth=2"
},

{
    "location": "man/Prerequisites/Interfaces.html#",
    "page": "Interfaces",
    "title": "Interfaces",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Prerequisites.Interfaces"
},

{
    "location": "man/Prerequisites/Interfaces.html#Interfaces-1",
    "page": "Interfaces",
    "title": "Interfaces",
    "category": "section",
    "text": "This submodule contains the generic functions that are extended by the package.Due to the multidispatch feature of Julia, generic functions can be extended by local methods for different types. However, a local definition of a method also claims a new generic function if the generic function is not imported to the current scope, thus ruins the definitions in other modules. Therefore, it is quite necessary to prefine the common generic functions in a separate module, so that other modules can extend them with their own by a simple import."
},

{
    "location": "man/Prerequisites/Interfaces.html#Hamiltonian.Prerequisites.Interfaces.:⊕",
    "page": "Interfaces",
    "title": "Hamiltonian.Prerequisites.Interfaces.:⊕",
    "category": "function",
    "text": "Direct sum.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Interfaces.html#Hamiltonian.Prerequisites.Interfaces.:⊗",
    "page": "Interfaces",
    "title": "Hamiltonian.Prerequisites.Interfaces.:⊗",
    "category": "function",
    "text": "Direct product.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Interfaces.html#Hamiltonian.Prerequisites.Interfaces.add!",
    "page": "Interfaces",
    "title": "Hamiltonian.Prerequisites.Interfaces.add!",
    "category": "function",
    "text": "Inplace addition.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Interfaces.html#Hamiltonian.Prerequisites.Interfaces.degree",
    "page": "Interfaces",
    "title": "Hamiltonian.Prerequisites.Interfaces.degree",
    "category": "function",
    "text": "Degree\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Interfaces.html#Hamiltonian.Prerequisites.Interfaces.dimension",
    "page": "Interfaces",
    "title": "Hamiltonian.Prerequisites.Interfaces.dimension",
    "category": "function",
    "text": "Dimension.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Interfaces.html#Hamiltonian.Prerequisites.Interfaces.dims",
    "page": "Interfaces",
    "title": "Hamiltonian.Prerequisites.Interfaces.dims",
    "category": "function",
    "text": "Dimensions.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Interfaces.html#Hamiltonian.Prerequisites.Interfaces.div!",
    "page": "Interfaces",
    "title": "Hamiltonian.Prerequisites.Interfaces.div!",
    "category": "function",
    "text": "Inplace division.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Interfaces.html#Hamiltonian.Prerequisites.Interfaces.expand",
    "page": "Interfaces",
    "title": "Hamiltonian.Prerequisites.Interfaces.expand",
    "category": "function",
    "text": "Get the expansion.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Interfaces.html#Hamiltonian.Prerequisites.Interfaces.index",
    "page": "Interfaces",
    "title": "Hamiltonian.Prerequisites.Interfaces.index",
    "category": "function",
    "text": "Index.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Interfaces.html#Hamiltonian.Prerequisites.Interfaces.inds",
    "page": "Interfaces",
    "title": "Hamiltonian.Prerequisites.Interfaces.inds",
    "category": "function",
    "text": "Indices.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Interfaces.html#Hamiltonian.Prerequisites.Interfaces.matrix",
    "page": "Interfaces",
    "title": "Hamiltonian.Prerequisites.Interfaces.matrix",
    "category": "function",
    "text": "Matrix representation.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Interfaces.html#Hamiltonian.Prerequisites.Interfaces.mul!",
    "page": "Interfaces",
    "title": "Hamiltonian.Prerequisites.Interfaces.mul!",
    "category": "function",
    "text": "Inplace multiplication.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Interfaces.html#Hamiltonian.Prerequisites.Interfaces.permute",
    "page": "Interfaces",
    "title": "Hamiltonian.Prerequisites.Interfaces.permute",
    "category": "function",
    "text": "Get the permutation.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Interfaces.html#Hamiltonian.Prerequisites.Interfaces.rank",
    "page": "Interfaces",
    "title": "Hamiltonian.Prerequisites.Interfaces.rank",
    "category": "function",
    "text": "Rank.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Interfaces.html#Hamiltonian.Prerequisites.Interfaces.sub!",
    "page": "Interfaces",
    "title": "Hamiltonian.Prerequisites.Interfaces.sub!",
    "category": "function",
    "text": "Inplace subtraction.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Interfaces.html#Hamiltonian.Prerequisites.Interfaces.vector",
    "page": "Interfaces",
    "title": "Hamiltonian.Prerequisites.Interfaces.vector",
    "category": "function",
    "text": "Vector representation.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Interfaces.html#Hamiltonian.Prerequisites.Interfaces.update!",
    "page": "Interfaces",
    "title": "Hamiltonian.Prerequisites.Interfaces.update!",
    "category": "function",
    "text": "Inplace Update.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Interfaces.html#Manual-1",
    "page": "Interfaces",
    "title": "Manual",
    "category": "section",
    "text": "Modules=[Interfaces]\nOrder=  [:module,:constant,:type,:macro,:function]"
},

{
    "location": "man/Prerequisites/TypeTraits.html#",
    "page": "Type traits",
    "title": "Type traits",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Prerequisites.TypeTraits"
},

{
    "location": "man/Prerequisites/TypeTraits.html#Type-traits-1",
    "page": "Type traits",
    "title": "Type traits",
    "category": "section",
    "text": "This module defines generic type traits that are useful to the package.Julia does not support multi-inheritance, which is sometimes not convenient. A way around this is to use traits, i.e. by utilizing the dispatch on certain (singleton) types known as traits to simulate multi-inheritance. Although this method cannot aviod small repetitive codes, it suits methods well that are complicated and lengthy."
},

{
    "location": "man/Prerequisites/TypeTraits.html#EfficientOperations-1",
    "page": "Type traits",
    "title": "EfficientOperations",
    "category": "section",
    "text": "EfficientOperations defines efficient operations such as ==/isequal, </isless, replace, etc, that ensure type stability.Type stability is the key of Julia to improve the code efficiency. However, it cannot be ensured in some unexpected cases, especially where an iterator is involved. For example, the following codes appears type unstable:function Base.:(==)(o1::AbstractType,o2::AbstractType)\n    n1,n2=o1|>typeof|>fieldcount,o2|>typeof|>fieldcount\n    n1==n2 ? all(getfield(o1,i)==getfield(o2,i) for i=1:n1) : false\nendMethods like above are common when we design abstract types, but they are not type stable. To get rid of it, the generated function trick can be used:@generated function Base.:(==)(o1::AbstractType,o2::AbstractType)\n    n1,n2=o1|>typeof|>fieldcount,o2|>typeof|>fieldcount\n    if n1==n2\n        expr=:(getfield(o1,1)==getfield(o2,1))\n        for i=2:fcount\n            expr=Expr(:&&,expr,:(getfield(o1,$i)==getfield(o2,$i)))\n        end\n        return expr\n    else\n        return :(false)\n    end\nendThen type stability can be ensured. We use this trick to implement the methods such as ==/isequal, </isless, replace, etc, with the trait efficientoperations::EfficientOperations. Other types can resort to these methods by passing efficientoperations as the first argument."
},

{
    "location": "man/Prerequisites/TypeTraits.html#MemoryOrder-1",
    "page": "Type traits",
    "title": "MemoryOrder",
    "category": "section",
    "text": "MemoryOrder provides the convertions, subtoind and indtosub, between a Cartesian index represented by a tuple and a linear index represented by an integer. C/C++ order or Fortran order can be specified, though the constant instances corder or forder of singleton types COrder and FOrder, which are both subtypes of the abstract type MemoryOrder."
},

{
    "location": "man/Prerequisites/TypeTraits.html#Hamiltonian.Prerequisites.TypeTraits.corder",
    "page": "Type traits",
    "title": "Hamiltonian.Prerequisites.TypeTraits.corder",
    "category": "constant",
    "text": "corder\n\nIndicate that the convertion between Cartesian index and linear index is using the C/C++ order.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/TypeTraits.html#Hamiltonian.Prerequisites.TypeTraits.efficientoperations",
    "page": "Type traits",
    "title": "Hamiltonian.Prerequisites.TypeTraits.efficientoperations",
    "category": "constant",
    "text": "efficientoperations\n\nIndicate that the efficient operations, i.e. \"==\"/\"isequal\", \"<\"/\"isless\" or \"replace\", will be used.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/TypeTraits.html#Hamiltonian.Prerequisites.TypeTraits.forder",
    "page": "Type traits",
    "title": "Hamiltonian.Prerequisites.TypeTraits.forder",
    "category": "constant",
    "text": "forder\n\nIndicate that the convertion between Cartesian index and linear index is using the Fortran order.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/TypeTraits.html#Hamiltonian.Prerequisites.TypeTraits.indtosub-Tuple{Tuple,Int64,Hamiltonian.Prerequisites.TypeTraits.FOrder}",
    "page": "Type traits",
    "title": "Hamiltonian.Prerequisites.TypeTraits.indtosub",
    "category": "method",
    "text": "indtosub(dims::Tuple,ind::Int,order::FOrder) -> Tuple\nindtosub(dims::Tuple,ind::Int,order::COrder) -> Tuple\n\nConvert an linear index to Cartesian index. Fortran-order or C-order can be assigned.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/TypeTraits.html#Hamiltonian.Prerequisites.TypeTraits.subtoind-Union{Tuple{N}, Tuple{Tuple{Vararg{Int64,N}},Tuple{Vararg{Int64,N}},FOrder}} where N",
    "page": "Type traits",
    "title": "Hamiltonian.Prerequisites.TypeTraits.subtoind",
    "category": "method",
    "text": "subtoind(dims::NTuple{N,Int},inds::NTuple{N,Int},order::FOrder) where N -> Int\nsubtoind(dims::NTuple{N,Int},inds::NTuple{N,Int},order::COrder) where N -> Int\n\nConvert an Cartesian index to linear index. Fortran-order or C-order can be assigned.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/TypeTraits.html#Base.:<-Tuple{Hamiltonian.Prerequisites.TypeTraits.EfficientOperations,Any,Any}",
    "page": "Type traits",
    "title": "Base.:<",
    "category": "method",
    "text": "<(::EfficientOperations,o1,o2) -> Bool\nisless(::EfficientOperations,o1,o2) -> Bool\n\nCompare two objects and judge whether the first is less than the second.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/TypeTraits.html#Base.:==-Tuple{Hamiltonian.Prerequisites.TypeTraits.EfficientOperations,Any,Any}",
    "page": "Type traits",
    "title": "Base.:==",
    "category": "method",
    "text": "==(::EfficientOperations,o1,o2) -> Bool\nisequal(::EfficientOperations,o1,o2) -> Bool\n\nCompare two objects and judge whether they are eqaul to each other.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/TypeTraits.html#Base.replace-Tuple{Hamiltonian.Prerequisites.TypeTraits.EfficientOperations,Any}",
    "page": "Type traits",
    "title": "Base.replace",
    "category": "method",
    "text": "replace(::EfficientOperations,o;kwargs...) -> typeof(o)\n\nReturn a copy of the input object with some of the field values replaced by the keyword arguments.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/TypeTraits.html#Manual-1",
    "page": "Type traits",
    "title": "Manual",
    "category": "section",
    "text": "Modules=[TypeTraits]\nOrder=  [:module,:constant,:type,:macro,:function]"
},

{
    "location": "man/Prerequisites/Factories.html#",
    "page": "Factories",
    "title": "Factories",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Prerequisites.Factoriespush!(LOAD_PATH,\"../../../../src/\")\nusing Hamiltonian.Prerequisites.Factories"
},

{
    "location": "man/Prerequisites/Factories.html#Factories-1",
    "page": "Factories",
    "title": "Factories",
    "category": "section",
    "text": "The aim of Factories is to provide tools to hack into Julia codes without knowing the details of their abstract syntax trees and regularize the mechanism to \"escape\" variables in Expr expressions, so that users can manipulate the existing codes, modify them and generate new ones in macros. In particular, a factory in this module means the representation of certain blocks of Julia codes by a usual Julia struct. This representation is much easier to comprehend than the canonical Expr representation and makes it far more convenient to define macros. In general, we propose the following requirements that any factory must satisfy:DECOMPOSITION - An Expr expression can be decomposed into its corresponding factory by the factory\'s constructor.\nCOMPOSITION - A factory can compose its corresponding Expr expression by calling itself.\nESCAPE - A variable should be or not be escaped in the composed Expr expression by a factory depends on predefined escape mechanisms.These three requirements also define the basic interfaces to interact with factories. In practice, we combine the second and third in a single interface, i.e. by passing an instance of certain concrete EscapeMechanism as the only argument of calling a factory, the needed Expr expression with variables correctly escaped can be obtained."
},

{
    "location": "man/Prerequisites/Factories.html#Escape-mechanisms-1",
    "page": "Factories",
    "title": "Escape mechanisms",
    "category": "section",
    "text": "We adopt Julia structs to denote escape mechanisms so that we can utilize Julia\'s multidispatch to implement different mechanisms whereas keeping the same interface."
},

{
    "location": "man/Prerequisites/Factories.html#EscapeMechanism-1",
    "page": "Factories",
    "title": "EscapeMechanism",
    "category": "section",
    "text": "EscapeMechanism is the abstract type for all concrete escape mechanisms."
},

{
    "location": "man/Prerequisites/Factories.html#Escaped-1",
    "page": "Factories",
    "title": "Escaped",
    "category": "section",
    "text": "Escaped has only one attribute:names::NTuple{N,Symbol} where N: the names of variables to be escapedApprently, a variable should be escaped if its name is in the names of an Escaped. This mechanism suits a factory whose variables should be unescaped by default."
},

{
    "location": "man/Prerequisites/Factories.html#UnEscaped-1",
    "page": "Factories",
    "title": "UnEscaped",
    "category": "section",
    "text": "UnEscaped also has only on attribute:names::NTuple{N,Symbol} where N: the names of variables not to be escapedObviously, on the contrary to Escaped, a variable should be escaped if its name is not in the names of an UnEscaped. This mechanism suits a factory whose variables should be escaped by default."
},

{
    "location": "man/Prerequisites/Factories.html#MixEscaped-1",
    "page": "Factories",
    "title": "MixEscaped",
    "category": "section",
    "text": "MixEscaped has two attributes:escaped::Escaped: the escaped part of the mixed mechanism\nunescaped::UnEscaped: the UnEscaped part of the mixed mechanismThis mechanism suits complex factories that parts of it suit the \"escaped\" mechanism while others suit the \"unescaped\" mechanism."
},

{
    "location": "man/Prerequisites/Factories.html#RawExpr-1",
    "page": "Factories",
    "title": "RawExpr",
    "category": "section",
    "text": "RawExpr has no attributes and it means \"raw expression without any variable escaped\". This mechanism is used for the print of all factories by default."
},

{
    "location": "man/Prerequisites/Factories.html#Concrete-factories-1",
    "page": "Factories",
    "title": "Concrete factories",
    "category": "section",
    "text": "Out of practical purposes, we implemente 7 kinds of factories, i.e. Inference, Argument, Parameter, Field, Block, FunctionFactory and TypeFactory, which represent a type inference, a function argument, a method or type parameter, a struct field, a begin ... end block, a function itself and a struct itself, respectively. Some of the basic methods making the above three requirements fulfilled with these types are based on the powerful functions defined in MacroTools.We want to give a remark that although the types and functions provided in this module helps a lot for the definition of macros, macros should not be abused. On the one hand, some macros may change the language specifications, which makes it hard to understand the codes, and even splits the community; on the one hand, macros usually increases the precompiling/jit time, which means enormous uses of macros in a module may lead to an extremely long load time. Besides, due to the limited ability of the author, the codes in this module are not optimal, which adds to the jit overhead. Any promotion that keeps the interfaces unchanged is welcomed."
},

{
    "location": "man/Prerequisites/Factories.html#Inference-1",
    "page": "Factories",
    "title": "Inference",
    "category": "section",
    "text": "An Inference has 3 attributes:head::Union{Symbol,Nothing}: the head of the type inference, which must be one of (nothing,:(::),:(<:),:curly)\nname::Union{Symbol,Nothing}: the name of the type inference\nparams::Union{Inference,Vector{Inference},Nothing}: the parameters of the type inferenceAll valid expressions representing type inferences can be passed to the constructor:Inference(:T)\nInference(:(<:Number))\nInference(:(Vector{T}))\nInference(:(Vector{Tuple{String,Int}}))\nInference(:(Type{<:Number}))On the other hand, you can use the macro @inference to construct an Inference directly from a type inference:@inference Vector{Tuple{String,Int}}note: Note\nInference is a recursive struct, i.e. it recursively decomposes a type inference until the final type inference is just a Symbol.\nWhen the input expression is a Symbol, the head and params attributes of the resulting Inference is nothing. Otherwise, its head is the same with that of the input expression, and the args of the input expression will be further decomposed, whose result will be stored in params.\nWhen the head of the input expression is :(<:), the params is an Inference whereas when the head of the input expression is :curly, the params is a Vector{Inference}.Inference uses the UnEscaped mechanism to escape variables, e.g.Inference(:(Vector{T}))(UnEscaped()) |> println\nInference(:(Vector{T}))(UnEscaped(:T)) |> println\nInference(:(Vector{T}))(UnEscaped(:Vector,:T)) |> println"
},

{
    "location": "man/Prerequisites/Factories.html#Argument-1",
    "page": "Factories",
    "title": "Argument",
    "category": "section",
    "text": "An Argument has 4 attributes:name::Union{Symbol,Nothing}: the name of the argument\ntype::Inference: the type inference of the argument\nslurp::Bool: whether the argument should be expanded by ...\ndefault::Any: the default value of the argument, nothing for those with no default valuesAll valid expressions representing the arguments of functions can be passed to the constructor:Argument(:arg)\nArgument(:(arg::ArgType))\nArgument(:(arg::ArgType...))\nArgument(:(arg::ArgType=default))Or you can use the macro @argument for a direct construction from an argument declaration:@argument arg::ArgType=defaultThe construction from such expressions is based on the the MacroTools.splitarg function.Argument uses the MixEscaped mechanism to escape variables, with the UnEscaped mechanism for type and Escaped mechanism for default, e.g.Argument(:(arg::Real=zero(Int)))(MixEscaped(UnEscaped(),Escaped(:zero,:Int))) |> printlnIt can be seen that the name of an argument will never be escaped, which is obvious since the name of a function argument is always local. By the way, the composition of an Argument expression is based on the MacroTools.combinearg function."
},

{
    "location": "man/Prerequisites/Factories.html#Parameter-1",
    "page": "Factories",
    "title": "Parameter",
    "category": "section",
    "text": "A Parameter has 2 attributes:name::Union{Symbol,Nothing}: the name of the parameter\ntype::Union{Inference,Nothing}: the type inference of the parameterAll expressions that represent type parameters or method parameters are allowed to be passed to the constructor:Parameter(:T)\nParameter(:(<:Number))\nParameter(:(T<:Number))The macro @parameter completes the construction directly from a parameter declaration:@parameter T<:Numbernote: Note\nWe use nothing to denote a missing name or type.\nTwo subtle situations of type/method parameters, e.g. MyType{T} and MyType{Int}, should be distinguished by Parameter(:T) and Parameter(:(<:Int)). In other words, MyType{Int} is in fact not supported. Indeed, Parameter(:Int) will treat :Int as the parameter name but not the parameter type.Parameter uses the UnEscaped mechanism to escape variables, too, e.g.Parameter(:(N<:Vector{T}))(UnEscaped(:T)) |> printlnAs is similar to Argument, the name of a method/type parameter will never be escaped because of its local scope."
},

{
    "location": "man/Prerequisites/Factories.html#Field-1",
    "page": "Factories",
    "title": "Field",
    "category": "section",
    "text": "A Field has 2 attributes:name::Symbol: the name of the field\ntype::Inference: the type inference of the fieldLegal expressions can be used to construct a Field instance by its constructor:Field(:field)\nField(:(field::FieldType))\nField(:(field::ParametricType{T}))The macro @field is also provided to help the construction directly from a field declaration:@field field::FieldTypeThe construction from these expressions is based on the MacroTools.splitarg function.Field uses the UnEscaped mechanism to escape variables as well, e.g.Field(:(field::Dict{N,D}))(UnEscaped(:N,:D)) |> printlnThe name of a struct will never be escaped either because it is a local variable tightly binding to a struct. It is noted that the composition of field expressions is based on the MacroTools.combinefield function."
},

{
    "location": "man/Prerequisites/Factories.html#Block-1",
    "page": "Factories",
    "title": "Block",
    "category": "section",
    "text": "A Block has only one attribute:body::Vector{Any}: the body of the begin ... end blockAny expression can be passed to the constructor of Block:Block(:(x=1))\nBlock(:(x=1;y=2))\nBlock(:(begin x=1 end))\nBlock(quote\n        x=1\n        y=2\n    end)Or you can construct a Block instance directly from any code by the macro @block:@block x=1 y=2The body of a block can also be extended by the push! function or the @push! macro.note: Note\nThe body of a Block is somewhat \"flattened\", i.e. it contains no begin ... end blocks. During the initialization, any such input block will be unblocked and added to the body part by part. So is the push! and @push! procedures.\nAll LineNumberNodes generated by the input codes will also be included in the block\'s body. However, you can use rmlines! or @rmlines! to remove them from the body of an existing Block, or use rmlines or @rmlines to get a copy with them removed in the body.Block uses the Escaped mechanism to escape variables. This is because variables in a block are often local ones and should not be escaped. Therefore, only those defined in other modules should be noted and escaped, which usually constitute the minority. For example,Block(:(x=1;y=2;z=Int[1,2,3]))(Escaped(:Int)) |> println"
},

{
    "location": "man/Prerequisites/Factories.html#FunctionFactory-1",
    "page": "Factories",
    "title": "FunctionFactory",
    "category": "section",
    "text": "A FunctionFactory has 7 attributes:name::Union{Symbol,Expr}: the name of the function\nparams::Vector{Inference}: the method parameters of the function\nargs::Vector{Argument}: the positional arguments of the function\nkwargs::Vector{Argument}: the keyword arguments of the function\nrtype::Inference: the return type of the function\nwhereparams::Vector{Parameter}: the method parameters specified by the where keyword\nbody::Block: the body of the functionAll expressions that represent functions are allowed to be passed to the constructor:FunctionFactory(:(f()=nothing))\nFunctionFactory(:(f(x)=x))\nFunctionFactory(:(f(x::Int,y::Int;choice::Function=sum)=choice(x,y)))\nFunctionFactory(:(f(x::T,y::T;choice::Function=sum) where T<:Number=choice(x,y)))\nFunctionFactory(:((f(x::T,y::T;choice::Function=sum)::T) where T<:Number=choice(x,y)))\nFunctionFactory(:(\n    function (f(x::T,y::T;choice::Function=sum)::T) where T<:Number\n        choice(x,y)\n    end\n))\nFunctionFactory(\n    quote\n        function (f(x::T,y::T;choice::Function=sum)::T) where T<:Number\n            choice(x,y)\n        end\n    end\n)Similarly, an instance can also be constructed from the macro @functionfactory:@functionfactory (f(x::T,y::T;choice::Function=sum)::T) where T<:Number=choice(x,y)The construction from such expressions are based on the MacroTools.splitdef function.note: Note\nSince Julia 0.7, the form MyType{D}(data::D) where D only appears in struct constructors, therefore, the attribute :params of a function factory is nonempty only when this factory aims to represent a struct constructor.\nUsually, the name of a function factory is a Symbol. However, if the factory aims to extend some methods of a function defined in another module, e.g., Base.eltype, the name will be an Expr.FunctionFactory adopts the MixEscaped mechanism to escape variables, with UnEscaped for params, args, kwargs, rtype and whereparams while Escaped for name and body. It is worth to emphasize that the name of a function factory belongs to the Escaped part. Therefore, when it is an Expr, it will never be escaped because an Expr cannot be a element of a NTuple{N,Symbol} where N. See examples,FunctionFactory(:(\n    (f(x::T,y::T;choice::Function=sum)::T) where T<:Number=max(x,y,choice(x,y))\n    ))(MixEscaped(UnEscaped(:T),Escaped(:f,:max,))) |> println\nFunctionFactory(:(\n    (f(x::T,y::T;choice::Function=sum)::T) where T<:Number=max(x,y,choice(x,y))\n    ))(MixEscaped(UnEscaped(:T),Escaped(:max))) |> printlnThe compositions of function expressions are based on the MacroTools.combinedef function.Other features include:Positional arguments can be added by addargs! or @addargs!\nKeyword arguments can be added by addkwargs! or @addkwargs!\nWhere parameters can be added by addwhereparams! or @addwhereparams!\nBody can be extended by extendbody! or @extendbody!"
},

{
    "location": "man/Prerequisites/Factories.html#TypeFactory-1",
    "page": "Factories",
    "title": "TypeFactory",
    "category": "section",
    "text": "A TypeFactory has 6 attributes:name::Symbol: the name of the struct\nmutable::Bool: whether or not the struct is mutable\nparams::Vector{Parameter}: the type parameters of the struct\nsupertype::Inference: the supertype of the struct\nfields::Vector{Field}: the fields of the struct\nconstructors::Vector{FunctionFactory}: the inner constructors of the structAny expression representing valid struct definitions can be passed to the constructor:TypeFactory(:(struct StructName end))\nTypeFactory(:(struct StructName{T} end))\nTypeFactory(:(struct Child{T} <: Parent{T} end))\nTypeFactory(:(\n    struct Child{T<:Number} <: Parent{T}\n        field1::T\n        field2::T\n    end\n))\nTypeFactory(\n    quote\n        struct Child{T<:Number} <: Parent{T}\n            field1::T\n            field2::T\n        end\n    end\n)Also, the macro @typefactory supports the construction directly from a type definition:@typefactory struct Child{T<:Number} <: Parent{T}\n                field1::T\n                field2::T\n                Child(field1::T,field2::T=zero(T)) where T=new{T}(field1,field2)\n            endThe construction from these expressions is based on the MacroTools.splitstructdef function.TypeFactory also uses the MixEscaped mechanism to escape variables, with the UnEscaped part for params, supertype and fields, the Escaped part for name, and both for constructors. For example,@typefactory(struct Child{T<:Number} <: Parent{T}\n    field::T\n    Child(field::T) where T=new{T}(field)\nend)(MixEscaped(UnEscaped(:T),Escaped(:Child))) |>printlnThe composition of a type expression is based on the MacroTools.combinestructdef function.Other features include:Fields can be added by addfields! or @addfields!\nType parameters can be added by addparams! or @addparams!\nInner constructors can be added by addconstructors! or @addconstructors!"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.FExpr",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.FExpr",
    "category": "constant",
    "text": "Factory expression types, which is defined as Union{Symbol,Expr}.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.rawexpr",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.rawexpr",
    "category": "constant",
    "text": "rawexpr\n\nIndicate that no variable in a factory should be escaped.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.Argument",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.Argument",
    "category": "type",
    "text": "Argument(name::Union{Symbol,Nothing},type::Inference,slurp::Bool,default::Any)\nArgument(;name::Union{Symbol,Nothing}=nothing,type::Inference=Inference(:Any),slurp::Bool=false,default::Any=nothing)\nArgument(expr::FExpr)\n\nThe struct to describe a argument of a function.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.Argument-Tuple{Hamiltonian.Prerequisites.Factories.RawExpr}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.Argument",
    "category": "method",
    "text": "(a::Argument)(em::RawExpr) -> Expr\n(a::Argument)(em::MixEscaped) -> Expr\n\nConvert an Argument to the Expr representation of the argument it describes.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.Block",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.Block",
    "category": "type",
    "text": "Block(parts::FExpr...)\n\nThe struct to describe a begin ... end block.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.Block-Tuple{Hamiltonian.Prerequisites.Factories.MixEscaped}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.Block",
    "category": "method",
    "text": "(b::Block)(em::RawExpr) -> Expr\n(b::Block)(em::Escaped) -> Expr\n(b::Block)(em::MixEscaped) -> Expr\n\nConvert a Block to the Expr representation of the begin ... end block it describes.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.Escaped",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.Escaped",
    "category": "type",
    "text": "Escaped(names::Symbol...)\n\nIndicate that symbols of a factory should be escaped if they are in names.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.Field",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.Field",
    "category": "type",
    "text": "Field(name::Symbol,type::Inference)\nField(;name::Symbol,type::FExpr=Inference(:Any))\nField(expr::FExpr)\n\nThe struct to describe a field of a struct.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.Field-Tuple{Union{RawExpr, #s14, #s13} where #s13<:Hamiltonian.Prerequisites.Factories.MixEscaped where #s14<:Hamiltonian.Prerequisites.Factories.UnEscaped}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.Field",
    "category": "method",
    "text": "(f::Field)(em::RawExpr) -> Expr\n(f::Field)(em::UnEscaped) -> Expr\n(f::Field)(em::MixEscaped) -> Expr\n\nConvert a Field to the Expr representation of the field it describes.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.FunctionFactory",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.FunctionFactory",
    "category": "type",
    "text": "FunctionFactory(name::FExpr,params::Vector{Inference},args::Vector{Argument},kwargs::Vector{Argument},rtype::Inference,whereparams::Vector{Parameter},body::Block)\nFunctionFactory(    ;name::FExpr,\n                    params::Vector{Inference}=Inference[],\n                    args::Vector{Argument}=Argument[],\n                    kwargs::Vector{Argument}=Argument[],\n                    rtype::Inference=Inference(:Any),\n                    whereparams::Vector{Parameter}=Parameter[],\n                    body::Block=Block()\n                    )\nFunctionFactory(expr::Expr)\n\nThe struct to describe a function.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.FunctionFactory-Tuple{Union{RawExpr, #s68} where #s68<:Hamiltonian.Prerequisites.Factories.MixEscaped}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.FunctionFactory",
    "category": "method",
    "text": "(ff::FunctionFactory)(em::RawExpr) -> Expr\n(ff::FunctionFactory)(em::MixEscaped) -> Expr\n\nConvert a FunctionFactory to the Expr representation of the function it describes.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.Inference",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.Inference",
    "category": "type",
    "text": "Inference(head::Union{Symbol,Nothing},name::Union{Symbol,Nothing},params::Union{Inference,Vector{Inference},Nothing})\nInference(;\n        head::Union{Symbol,Nothing}=nothing,\n        name::Union{Symbol,Nothing}=nothing,\n        params::Union{Inference,Vector{Inference},Nothing}=nothing,\n        )\nInference(expr::FExpr)\n\nThe struct to describe a type inference.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.Inference-Tuple{Hamiltonian.Prerequisites.Factories.MixEscaped}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.Inference",
    "category": "method",
    "text": "(i::Inference)(em::RawExpr) -> FExpr\n(i::Inference)(em::UnEscaped) -> FExpr\n(i::Inference)(em::MixEscaped) -> FExpr\n\nConvert a Inference to the Expr representation of the type inference it describes.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.MixEscaped",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.MixEscaped",
    "category": "type",
    "text": "MixEscaped(escaped::Escaped)\nMixEscaped(unescaped::UnEscaped)\nMixEscaped(escaped::Escaped,unescaped::UnEscaped)\nMixEscaped(unescaped::UnEscaped,escaped::Escaped)\n\nIndicate that some parts of a factory use the Escaped mechanism while other parts use the UnEscaped mechanism.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.Parameter",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.Parameter",
    "category": "type",
    "text": "Parameter(name::Union{Symbol,Nothing},type::Union{Inference,Nothing})\nParameter(;name::Union{Symbol,Nothing}=nothing,type::Union{Inference,Nothing}=nothing)\nParameter(expr::FExpr)\n\nThe struct to describe a parameter of a function or a type.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.Parameter-Tuple{Union{RawExpr, #s14, #s13} where #s13<:Hamiltonian.Prerequisites.Factories.MixEscaped where #s14<:Hamiltonian.Prerequisites.Factories.UnEscaped}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.Parameter",
    "category": "method",
    "text": "(p::Parameter)(em::RawExpr) -> FExpr\n(p::Parameter)(em::UnEscaped) -> FExpr\n(p::Parameter)(em::MixEscaped) -> FExpr\n\nConvert a Parameter to the Expr representation of the parameter it describes.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.RawExpr",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.RawExpr",
    "category": "type",
    "text": "Raw expression without any variable escaped.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.TypeFactory",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.TypeFactory",
    "category": "type",
    "text": "TypeFactory(name::Symbol,mutable::Bool,params::Vector{Parameter},supertype::Inference,fields::Vector{Field},constructors::Vector{FunctionFactory})\nTypeFactory(    ;name::Symbol,\n                mutable::Bool=false,\n                params::Vector{Parameter}=Parameter[],\n                supertype::Inference=Inference(:Any),\n                fields::Vector{Field}=Field[],\n                constructors::Vector{FunctionFactory}=FunctionFactory[],\n                )\nTypeFactory(expr::Expr)\n\nThe struct to describe a struct.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.TypeFactory-Tuple{Union{RawExpr, #s68} where #s68<:Hamiltonian.Prerequisites.Factories.MixEscaped}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.TypeFactory",
    "category": "method",
    "text": "(tf::TypeFactory)(em::RawExpr) -> Expr\n(tf::TypeFactory)(em::MixEscaped) -> Expr\n\nConvert a TypeFactory to the Expr representation of the struct it describes.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.UnEscaped",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.UnEscaped",
    "category": "type",
    "text": "UnEscaped(names::Symbol...)\n\nIIndicate that symbols of a factory should be escaped if they are not in names.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.@addargs!-Tuple{Any,Vararg{Union{Expr, Symbol},N} where N}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.@addargs!",
    "category": "macro",
    "text": "@addargs! ff args::FExpr...\n\nAdd a couple of positional arguments to a function factory.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.@addconstructors!-Tuple{Any,Vararg{Expr,N} where N}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.@addconstructors!",
    "category": "macro",
    "text": "@addconstructors! tf constructors::Expr...\n\nAdd a couple of constructors to a type factory.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.@addfields!-Tuple{Any,Vararg{Union{Expr, Symbol},N} where N}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.@addfields!",
    "category": "macro",
    "text": "@addfields! tf fields::FExpr...\n\nAdd a couple of fields to a type factory.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.@addkwargs!-Tuple{Any,Vararg{Union{Expr, Symbol},N} where N}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.@addkwargs!",
    "category": "macro",
    "text": "@addkwargs! ff kwargs::FExpr...\n\nAdd a couple of keyword arguments to a function factory.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.@addparams!-Tuple{Any,Vararg{Union{Expr, Symbol},N} where N}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.@addparams!",
    "category": "macro",
    "text": "@addparams! f params::FExpr...\n\nAdd a couple of method parameters to a function factory or a type factory.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.@addwhereparams!-Tuple{Any,Vararg{Union{Expr, Symbol},N} where N}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.@addwhereparams!",
    "category": "macro",
    "text": "@addwhereparams! ff whereparams::FExpr...\n\nAdd a couple of method parameters to a function factory or a type factory.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.@argument-Tuple{Union{Expr, Symbol}}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.@argument",
    "category": "macro",
    "text": "@argument expr::FExpr\n\nConstruct an Argument directly from an argument statement.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.@block-Tuple{Vararg{Union{Expr, Symbol},N} where N}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.@block",
    "category": "macro",
    "text": "@block parts::FExpr...\n\nConstruct a Block directly from a begin ... end block definition.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.@extendbody!-Tuple{Any,Vararg{Union{Expr, Symbol},N} where N}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.@extendbody!",
    "category": "macro",
    "text": "@extendbody! ff parts::FExpr...\n\nExtend the body of a function factory.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.@field-Tuple{Union{Expr, Symbol}}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.@field",
    "category": "macro",
    "text": "@field expr::FExpr\n\nConstruct a Field directly from a field statement.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.@functionfactory-Tuple{Expr}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.@functionfactory",
    "category": "macro",
    "text": "@functionfactory expr::FExpr\n\nConstruct a FunctionFactory directly from a function definition.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.@inference-Tuple{Union{Expr, Symbol}}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.@inference",
    "category": "macro",
    "text": "@inference expr::FExpr\n\nConstruct an Inference directly from a type inference.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.@parameter-Tuple{Union{Expr, Symbol}}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.@parameter",
    "category": "macro",
    "text": "@parameter expr::FExpr\n\nConstruct a Parameter directly from an parameter statement.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.@push!-Tuple{Any,Vararg{Union{Expr, Symbol},N} where N}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.@push!",
    "category": "macro",
    "text": "@push! b parts::FExpr...\n\nPush other parts into the body of a block.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.@rmlines!-Tuple{Expr}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.@rmlines!",
    "category": "macro",
    "text": "@rmlines! b::Expr\n\nRemove line number nodes in the body of a block.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.@rmlines-Tuple{Expr}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.@rmlines",
    "category": "macro",
    "text": "@rmlines b::Expr\n\nReturn a copy of a block with the line number nodes removed.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.@typefactory-Tuple{Expr}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.@typefactory",
    "category": "macro",
    "text": "@typefactory expr::Expr\n\nConstruct a TypeFactory directly from a type definition.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.addargs!-Tuple{Hamiltonian.Prerequisites.Factories.FunctionFactory}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.addargs!",
    "category": "method",
    "text": "addargs!(ff::FunctionFactory,args::Argument...) -> FunctionFactory\naddargs!(ff::FunctionFactory,args::FExpr...) -> FunctionFactory\n\nAdd a couple of positional arguments to a function factory.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.addconstructors!-Tuple{Hamiltonian.Prerequisites.Factories.TypeFactory}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.addconstructors!",
    "category": "method",
    "text": "addconstructors!(tf::TypeFactory,constructors::FunctionFactory...) -> TypeFactory\naddconstructors!(tf::TypeFactory,constructors::Expr...) -> TypeFactory\n\nAdd a couple of constructors to a type factory.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.addfields!-Tuple{Hamiltonian.Prerequisites.Factories.TypeFactory}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.addfields!",
    "category": "method",
    "text": "addfields!(tf::TypeFactory,fields::Field...) -> TypeFactory\naddfields!(tf::TypeFactory,fields::FExpr...) -> TypeFactory\n\nAdd a couple of fields to a type factory.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.addkwargs!-Tuple{Hamiltonian.Prerequisites.Factories.FunctionFactory}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.addkwargs!",
    "category": "method",
    "text": "addkwargs!(ff::FunctionFactory,kwargs::Argument...) -> FunctionFactory\naddkwargs!(ff::FunctionFactory,kwargs::FExpr...) -> FunctionFactory\n\nAdd a couple of keyword arguments to a function factory.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.addparams!-Tuple{Hamiltonian.Prerequisites.Factories.FunctionFactory}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.addparams!",
    "category": "method",
    "text": "addparams!(f::FunctionFactory,params::Inference...) ->FunctionFactory\naddparams!(f::FunctionFactory,params::FExpr...) -> FunctionFactory\naddparams!(f::TypeFactory,params::Parameter...) -> TypeFactory\naddparams!(f::TypeFactory,params::FExpr...) -> TypeFactory\n\nAdd a couple of parameters to a function factory or a type factory.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.addwhereparams!-Tuple{Hamiltonian.Prerequisites.Factories.FunctionFactory}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.addwhereparams!",
    "category": "method",
    "text": "addwhereparams!(ff::FunctionFactory,whereparams::Parameter...) -> FunctionFactory\naddwhereparams!(ff::FunctionFactory,whereparams::FExpr...) -> FunctionFactory\n\nAdd a couple of method where parameters to a function factory or a type factory.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.escape-Tuple{Any,Hamiltonian.Prerequisites.Factories.RawExpr}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.escape",
    "category": "method",
    "text": "escape(expr,::RawExpr) -> Any\nescape(expr,::Escaped) -> Any\nescape(expr,::UnEscaped) -> Any\nescape(expr::Symbol,em::Escaped) -> FExpr\nescape(expr::Expr,em::Escaped) -> Expr\nescape(expr::Symbol,em::UnEscaped) -> FExpr\nescape(expr::Expr,em::UnEscaped) -> Expr\n\nEscape the variables in the input expression.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.extendbody!-Tuple{Hamiltonian.Prerequisites.Factories.FunctionFactory}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.extendbody!",
    "category": "method",
    "text": "extendbody!(ff::FunctionFactory,parts::FExpr...) -> FunctionFactory\nextendbody!(ff::FunctionFactory,parts::Block...) -> FunctionFactory\n\nExtend the body of a function factory.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.rmlines!-Tuple{Hamiltonian.Prerequisites.Factories.Block}",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.rmlines!",
    "category": "method",
    "text": "rmlines!(b::Block) -> Block\n\nRemove line number nodes in the body of a block.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#MacroTools.rmlines-Tuple{Hamiltonian.Prerequisites.Factories.Block}",
    "page": "Factories",
    "title": "MacroTools.rmlines",
    "category": "method",
    "text": "rmlines(b::Block) -> Block\n\nReturn a copy of a block with the line number nodes removed.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.AbstractFactory",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.AbstractFactory",
    "category": "type",
    "text": "AbstractFactory\n\nAbstract type for all concrete factories.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Hamiltonian.Prerequisites.Factories.EscapeMechanism",
    "page": "Factories",
    "title": "Hamiltonian.Prerequisites.Factories.EscapeMechanism",
    "category": "type",
    "text": "Abstract escape mechanism.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Base.:==-Union{Tuple{F}, Tuple{F,F}} where F<:Hamiltonian.Prerequisites.Factories.AbstractFactory",
    "page": "Factories",
    "title": "Base.:==",
    "category": "method",
    "text": "==(f1::F,f2::F) where F<:AbstractFactory -> Bool\nisequal(f1::F,f2::F) where F<:AbstractFactory -> Bool\n\nOverloaded equivalent operator.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Base.push!-Tuple{Hamiltonian.Prerequisites.Factories.Block}",
    "page": "Factories",
    "title": "Base.push!",
    "category": "method",
    "text": "push!(b::Block,parts::FExpr...) -> Block\npush!(b::Block,parts::Block...) -> Block\n\nPush other parts into the body of a block.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Base.replace-Tuple{Hamiltonian.Prerequisites.Factories.AbstractFactory}",
    "page": "Factories",
    "title": "Base.replace",
    "category": "method",
    "text": "replace(f::AbstractFactory;kwargs...) -> typeof(f)\n\nReturn a copy of a concrete AbstractFactory with some of the field values replaced by the keyword arguments.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Base.show-Tuple{IO,Hamiltonian.Prerequisites.Factories.AbstractFactory}",
    "page": "Factories",
    "title": "Base.show",
    "category": "method",
    "text": "show(io::IO,f::AbstractFactory)\n\nShow a concrete AbstractFactory.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Factories.html#Manual-1",
    "page": "Factories",
    "title": "Manual",
    "category": "section",
    "text": "Modules=[Factories]\nOrder=  [:module,:constant,:type,:macro,:function]"
},

{
    "location": "man/Prerequisites/CompositeStructures.html#",
    "page": "Composite structures",
    "title": "Composite structures",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Prerequisites.CompositeStructures"
},

{
    "location": "man/Prerequisites/CompositeStructures.html#Composite-structures-1",
    "page": "Composite structures",
    "title": "Composite structures",
    "category": "section",
    "text": "In principle, Julia is not an object-oriented programming language. For example, only abstract types can be inherited so that subtype cannot inherit fields from their parents. Therefore, Julia prefers composition over inheritance. However, to make a new concrete type behaves much alike another one, tedious reputitions of redifining the generic interfaces are usually not avoidable, especially for the basic types in Julia base. In this module, we implement three such composited types, CompositeTuple, CompositeVector and CompositeDict, for the sake of future usages."
},

{
    "location": "man/Prerequisites/CompositeStructures.html#CompositeTuple-and-CompositeNTuple-1",
    "page": "Composite structures",
    "title": "CompositeTuple and CompositeNTuple",
    "category": "section",
    "text": "A composite tuple (ntuple) is a tuple (ntuple) that is implemented by including an ordinary Tuple (NTuple) as one of its attributes with the name :contents.To take full advantages of the Julia base, the following interfaces are defined:inquiry of info: length, eltype, hash\ncomparison between objects: ==, isequal\nobtainment of old elements: getindex\niteration: iterate, keys, values, pairs\nconstruction of new objects: reverseComposite ntuples are suited for the situations where other attributes are not affected by the modification of the elements. Note that arithmatic operations and logical operations excluding == and isequal are not supported. Besides, a composite ntuple is not a tuple since Julia has no abstract tuples."
},

{
    "location": "man/Prerequisites/CompositeStructures.html#CompositeVector-1",
    "page": "Composite structures",
    "title": "CompositeVector",
    "category": "section",
    "text": "A composite vector is a vector that is implemented by including an ordinary Vector as one of its attributes with the name :contents.To take full advantages of the Julia base, the following interfaces are redefined:inquiry of info: size, length\ncomparison between objects: ==, isequal\nobtainment of old elements: getindex\nmodification of old elements: setindex!\naddition of new elements: push!, pushfirst!, insert!, append!, prepend!\nremoval of old elements: splice!, deleteat!, pop!, popfirst!, empty!\nconstruction of new objects: empty, reverse, similar\niteration: iterate, keys, values, pairsComposite vectors are suited for the situations where other attributes are not affected by the modification of the elements. Note that arithmatic operations and logical operations excluding == and isequal are not supported."
},

{
    "location": "man/Prerequisites/CompositeStructures.html#CompositeDict-1",
    "page": "Composite structures",
    "title": "CompositeDict",
    "category": "section",
    "text": "A composite dict is a dict that is implemented by including an ordinary Dict as one of its attributes with the name :contents.To take full advantages of the Julia base, the following interfaces are redefined:inquiry of info: isempty, length, haskey, in\ncomparison between objects: ==, isequal\nobtainment of old elements: get, getkey, getindex\nmodification and addition of elements: push!, get!, setindex!\nremoval of old elements: pop!, delete!, empty!\nconstruction of new objects: merge, empty\niteration: iterate, keys, values, pairsAs is similar to composite vectors, composite dicts are suited for the situations where other attributes are not affected by the modification of the elements."
},

{
    "location": "man/Prerequisites/CompositeStructures.html#Hamiltonian.Prerequisites.CompositeStructures.CompositeDict",
    "page": "Composite structures",
    "title": "Hamiltonian.Prerequisites.CompositeStructures.CompositeDict",
    "category": "type",
    "text": "CompositeDict{K,V}\n\nA composite dict is a dict that is implemented by including an ordinary Dict as one of its attributes with the name :contents.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/CompositeStructures.html#Hamiltonian.Prerequisites.CompositeStructures.CompositeNTuple",
    "page": "Composite structures",
    "title": "Hamiltonian.Prerequisites.CompositeStructures.CompositeNTuple",
    "category": "type",
    "text": "CompositeNTuple{N,T}\n\nA composite ntuple is a ntuple that is implemented by including an ordinary NTuple as one of its attributes with the name :contents.\n\nAlias for CompositeTuple{NTuple{N,T}}.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/CompositeStructures.html#Hamiltonian.Prerequisites.CompositeStructures.CompositeTuple",
    "page": "Composite structures",
    "title": "Hamiltonian.Prerequisites.CompositeStructures.CompositeTuple",
    "category": "type",
    "text": "CompositeTuple{T<:Tuple}\n\nA composite tuple is a tuple that is implemented by including an ordinary Tuple as one of its attributes with the name :contents.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/CompositeStructures.html#Hamiltonian.Prerequisites.CompositeStructures.CompositeVector",
    "page": "Composite structures",
    "title": "Hamiltonian.Prerequisites.CompositeStructures.CompositeVector",
    "category": "type",
    "text": "CompositeVector{T}\n\nA composite vector is a vector that is implemented by including an ordinary Vector as one of its attributes with the name :contents.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/CompositeStructures.html#Manul-1",
    "page": "Composite structures",
    "title": "Manul",
    "category": "section",
    "text": "Modules=[CompositeStructures]\nOrder=  [:module,:constant,:type,:macro,:function]"
},

{
    "location": "man/Prerequisites/Trees.html#",
    "page": "Trees",
    "title": "Trees",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Prerequisites.Treespush!(LOAD_PATH,\"../../../../src/\")\nusing Hamiltonian.Prerequisites.Trees"
},

{
    "location": "man/Prerequisites/Trees.html#Trees-1",
    "page": "Trees",
    "title": "Trees",
    "category": "section",
    "text": "The aim of this module is to represent the standard tree structure in efficiency-non-sensitive cases. Please note that the default implementation of tree methods are far from optimal in efficiency. Therefore, please DO NOT use it if you need an efficient tree for addition, deletion, sort and inquiry. This module of codes apply only when the structure of tree matters but not the efficiency."
},

{
    "location": "man/Prerequisites/Trees.html#AbstractTree-1",
    "page": "Trees",
    "title": "AbstractTree",
    "category": "section",
    "text": "AbstractTree{N,D} is the abstract type for all concrete trees. By design, it has two type parameters:N: the type of the tree\'s node\nD: the type of the tree\'s dataTo fully utilize the methods designed for a tree structure, in our protocol, a concrete subtype must implement the following methods:inquiry related methods\nroot(tree::AbstractTree{N,D}) where {N,D} -> Union{N,Nothing}\nGet a tree\'s root node (nothing for empty trees)\nhaskey(tree::AbstractTree{N,D},node::N) where {N,D} -> Bool\nCheck whether a node is in a tree.\nlength(tree::AbstractTree) -> Int\nGet the number of a tree\'s nodes.\nparent(tree::AbstractTree{N,D},\n       node::N,\n       superparent::Union{N,Nothing}=nothing\n       ) where {N,D} -> Union{N,Nothing}\nGet the parent of a tree\'s node or return superparent when the input node is the tree\'s root.\nchildren(tree::AbstractTree{N,D},node::N) where {N,D} -> Vector{N}\nGet the children of a tree\'s node.\nstructure modification related methods\naddnode!(tree::AbstractTree{N,D},\n         parent::Union{N,Nothing},\n         node::N\n         ) where {N,D}\nUpdate the structure of a tree by adding a node. When the parent is nothing, the input tree must be empty and the input node becomes the tree\'s root.\ndeletenode!(tree::AbstractTree{N,D},node::N) where {N,D}\nUpdate the structure of a tree by deleting a node.\nindex related methods\ngetindex(tree::AbstractTree{N,D},node::N) where {N,D} -> D\nGet the data of a tree\'s node\nsetindex!(tree::AbstractTree{N,D},node::N,data::D) where {N,D}\nSet the data of a tree\'s node.Based on these methods, we implement several generic functions for inquiries and manipulationsinquiry for type parameters: keytype, valtype, eltype\nexpansion over nodes/data-records: keys, values, pairs\ninquiry for info of nodes: isleaf, level\ninquiry for nodes: ancestor, descendants, siblings, leaves\nmodification: push!, append!, delete!, empty!And optionally, when a subtype implement the following method,empty(tree::AbstractTree) -> typeof(tree)which constructs an empty tree of the same type with the input one, two more more methods are supported:subtree: Get a subtree starting from a node.\nmove!: Move a subtree to a new position."
},

{
    "location": "man/Prerequisites/Trees.html#TreeCore-and-SimpleTree-1",
    "page": "Trees",
    "title": "TreeCore and SimpleTree",
    "category": "section",
    "text": "To implement all the prerequisites listed above costs a bit efforts. We provide two lazy ways to get over this:Inheritance AbstractTree with TREECORE::TreeCore as the last attribute\nInclusion an attribute which is an instance of SimpleTree"
},

{
    "location": "man/Prerequisites/Trees.html#TreeCore-1",
    "page": "Trees",
    "title": "TreeCore",
    "category": "section",
    "text": "TreeCore{N,D}, as the literal meaning indicates, is the core of a tree. It encapsulates all the data structures needed by the default implementation, which constains 4 attributes:root::N: the tree\'s root node\ncontents::Dict{N,D}: the tree\'s (node,data) pairs\nparent::Dict{N,N}: records of the parent of each of the tree\'s nodes\nchildren::Dict{N,Vector{N}}: records of the children of each of the tree\'s nodesAs above, the first lazy way is to include this struct with the special name :TREECORE in your concrete subtype as the last attribute. This process can be even lazier, in that we provide a macro @tree to decorate your \"raw\" struct automatically, e.g.@tree(struct SimpleSubTree end)\n@tree(struct SubTreeWithTreeParameters end,\n      {N<:AbstractString,D<:Number}\n      )\n@tree(struct SubTreeWithCertainTreeParameters end,\n      {<:String,<:Int}\n      )\n@tree(struct SubTreeWithFields info::Vector{Int} end,\n      {N<:AbstractString,D<:Number}\n      )\n@tree(struct SubTreeWithParametricFields{T} info::Vector{T} end,\n      {N<:AbstractString,D<:Number}\n      )\n@tree(struct SubTreeWithOverlappedParametricFields{N} info::Vector{N} end,\n      {N<:AbstractString,D<:Number}\n      )"
},

{
    "location": "man/Prerequisites/Trees.html#SimpleTree-1",
    "page": "Trees",
    "title": "SimpleTree",
    "category": "section",
    "text": "SimpleTree{N,D} is the minimum struct that implements all the default tree methods. You can include an instance of it as an attribute in your own type to utilize all the tree methods."
},

{
    "location": "man/Prerequisites/Trees.html#Hamiltonian.Prerequisites.Trees.treedepth",
    "page": "Trees",
    "title": "Hamiltonian.Prerequisites.Trees.treedepth",
    "category": "constant",
    "text": "treedepth\n\nIndicate that the iteration over a tree is depth-first.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Hamiltonian.Prerequisites.Trees.treewidth",
    "page": "Trees",
    "title": "Hamiltonian.Prerequisites.Trees.treewidth",
    "category": "constant",
    "text": "treewidth\n\nIndicate that the iteration over a tree is width-first.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Hamiltonian.Prerequisites.Trees.AbstractTree",
    "page": "Trees",
    "title": "Hamiltonian.Prerequisites.Trees.AbstractTree",
    "category": "type",
    "text": "AbstractTree{Node,Data}\n\nAbstract type for all concrete trees.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Hamiltonian.Prerequisites.Trees.SimpleTree",
    "page": "Trees",
    "title": "Hamiltonian.Prerequisites.Trees.SimpleTree",
    "category": "type",
    "text": "SimpleTree{N,D}() where {N,D}\n\nThe minimum tree structure that implements all the default tree methods.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Hamiltonian.Prerequisites.Trees.TreeCore",
    "page": "Trees",
    "title": "Hamiltonian.Prerequisites.Trees.TreeCore",
    "category": "type",
    "text": "TreeCore()\n\nThe core of a tree.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Hamiltonian.Prerequisites.Trees.@tree",
    "page": "Trees",
    "title": "Hamiltonian.Prerequisites.Trees.@tree",
    "category": "macro",
    "text": "@tree structdef treeparams::Union{Expr,Nothing}=nothing\n\nDecorate a \"raw\" struct to be a subtype of AbstractTree.\n\nnote: Note\nA \"raw\" struct means:\nIt has no explicit supertype;\nIt has no inner constructor;\nIt has no attribute :TREECORE.\nThe keytype and valtype can be assigned by the argument treeparams in the form {keytype,valtype}.\nWhen the formal argument names of keytype and valtype are not assigned, they can be automatically generated by the functioin gensym. For example, all of the structs after the decration by the following codes\n@tree(struct TreeWithWrongTypeParameterNames{N} info::Vector{N} end)\n@tree(struct TreeWithWrongTypeParameterNames{N} info::Vector{N} end,\n      {::String,::Int}\n      )\n@tree(struct TreeWithWrongTypeParameterNames{N} info::Vector{N} end,\n      {<:AbstractString,<:Number}\n      )\nwill have three type parameters.\nWhen the formal argument names of keytype and valtype overlap with those of the raw struct type parameters, the duplicates will be considered as the same. For example, the decorated struct SubTreeWithOverlappedParametricFields by the following code\n@tree(struct TreeWithOverlappedParametricFields{N} info::Vector{N} end,\n      {N<:AbstractString,D<:Number}\n      )\nonly has two type parameters N<:AbstractString and D<:Number, where the N in the info::Vector{N} is the same N with that in the decorated attribute TREECORE::TreeCore{N,D}.\nWhen the formal argument names of keytype and valtype have no intersection with those of the raw struct type parameters, the type parameters of the decorated struct will be just extended by keytype and valtype. For example, the decorated struct SubTreeWithParametricFields by the following code\n@tree(struct TreeWithParametricFields{T} info::Vector{T} end,\n      {N<:AbstractString,D<:Number}\n      )\nhave 3 type parameters, T, N<:AbstractString and D<:Number.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Hamiltonian.Prerequisites.Trees.addnode!-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},N}} where D where N",
    "page": "Trees",
    "title": "Hamiltonian.Prerequisites.Trees.addnode!",
    "category": "method",
    "text": "addnode!(tree::AbstractTree{N,D},node::N) where {N,D} -> typeof(tree)\naddnode!(tree::AbstractTree{N,D},::Nothing,node::N) where {N,D} -> typeof(tree)\naddnode!(tree::AbstractTree{N,D},parent::N,node::N) where {N,D} -> typeof(tree)\n\nUpdate the structure of a tree by adding a node. When the parent is nothing, the input tree must be empty and the input node becomes the tree\'s root.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Hamiltonian.Prerequisites.Trees.ancestor-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},N}, Tuple{AbstractTree{N,D},N,Int64}} where D where N",
    "page": "Trees",
    "title": "Hamiltonian.Prerequisites.Trees.ancestor",
    "category": "method",
    "text": "ancestor(tree::AbstractTree{N,D},node::N,generation::Int=1) where {N,D} -> N\n\nGet the ancestor of a tree\'s node of the n-th generation.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Hamiltonian.Prerequisites.Trees.children-Tuple{Hamiltonian.Prerequisites.Trees.AbstractTree}",
    "page": "Trees",
    "title": "Hamiltonian.Prerequisites.Trees.children",
    "category": "method",
    "text": "children(tree::AbstractTree) -> Vector{keytype(tree)}\nchildren(tree::AbstractTree,::Nothing) -> Vector{keytype(tree)}\nchildren(tree::AbstractTree{N,D},node::N) where {N,D} -> Vector{N}\n\nGet the children of a tree\'s node.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Hamiltonian.Prerequisites.Trees.deletenode!-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},N}} where D where N",
    "page": "Trees",
    "title": "Hamiltonian.Prerequisites.Trees.deletenode!",
    "category": "method",
    "text": "deletenode!(tree::AbstractTree{N,D},node::N) where {N,D} -> typeof(tree)\n\nUpdate the structure of a tree by deleting a node.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Hamiltonian.Prerequisites.Trees.descendants-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},N}, Tuple{AbstractTree{N,D},N,Int64}} where D where N",
    "page": "Trees",
    "title": "Hamiltonian.Prerequisites.Trees.descendants",
    "category": "method",
    "text": "descendants(tree::AbstractTree{N,D},node::N,generation::Int=1) where {N,D} -> Vector{N}\n\nGet the descendants of a tree\'s node of the n-th generation.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Hamiltonian.Prerequisites.Trees.isleaf-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},N}} where D where N",
    "page": "Trees",
    "title": "Hamiltonian.Prerequisites.Trees.isleaf",
    "category": "method",
    "text": "isleaf(tree::AbstractTree{N,D},node::N) where{N,D} -> Bool\n\nJudge whether a tree\'s node is a leaf (a node without children) or not.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Hamiltonian.Prerequisites.Trees.leaves-Tuple{Hamiltonian.Prerequisites.Trees.AbstractTree}",
    "page": "Trees",
    "title": "Hamiltonian.Prerequisites.Trees.leaves",
    "category": "method",
    "text": "leaves(tree::AbstractTree) -> Vector{keytype(tree)}\n\nGet a tree\'s leaves.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Hamiltonian.Prerequisites.Trees.level-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},N}} where D where N",
    "page": "Trees",
    "title": "Hamiltonian.Prerequisites.Trees.level",
    "category": "method",
    "text": "level(tree::AbstractTree{N,D},node::N) where {N,D} -> Int\n\nGet the level of tree\'s node.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Hamiltonian.Prerequisites.Trees.move!-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},N,N}} where D where N",
    "page": "Trees",
    "title": "Hamiltonian.Prerequisites.Trees.move!",
    "category": "method",
    "text": "move!(tree::AbstractTree{N,D},node::N,parent::N) where {N,D} -> typeof(tree)\n\nMove a subtree to a new position.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Hamiltonian.Prerequisites.Trees.parent-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},N}, Tuple{AbstractTree{N,D},N,Union{Nothing, N}}} where D where N",
    "page": "Trees",
    "title": "Hamiltonian.Prerequisites.Trees.parent",
    "category": "method",
    "text": "parent(tree::AbstractTree{N,D},node::N,superparent::Union{N,Nothing}=nothing) where {N,D} -> Union{N,Nothing}\n\nGet the parent of a tree\'s node. When node is the tree\'s root, return superparent.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Hamiltonian.Prerequisites.Trees.root-Tuple{Hamiltonian.Prerequisites.Trees.AbstractTree}",
    "page": "Trees",
    "title": "Hamiltonian.Prerequisites.Trees.root",
    "category": "method",
    "text": "root(tree::AbstractTree) -> Union{keytype(tree),Nothing}\n\nGet a tree\'s root node.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Hamiltonian.Prerequisites.Trees.siblings-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},N}} where D where N",
    "page": "Trees",
    "title": "Hamiltonian.Prerequisites.Trees.siblings",
    "category": "method",
    "text": "siblings(tree::AbstractTree{N,D},node::N) where{N,D} -> Vector{N}\n\nGet the siblings (other nodes sharing the same parent) of a tree\'s node.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Hamiltonian.Prerequisites.Trees.subtree-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},N}} where D where N",
    "page": "Trees",
    "title": "Hamiltonian.Prerequisites.Trees.subtree",
    "category": "method",
    "text": "subtree(tree::AbstractTree{N,D},node::N) where{N,D} -> typeof(tree)\n\nGet a subtree whose root is node.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Base.:==-Union{Tuple{TC}, Tuple{TC,TC}} where TC<:Hamiltonian.Prerequisites.Trees.TreeCore",
    "page": "Trees",
    "title": "Base.:==",
    "category": "method",
    "text": "==(tc1::TC,tc2::TC) where TC<:TreeCore -> Bool\nisequal(tc1::TC,tc2::TC) where TC<:TreeCore -> Bool\n\nOverloaded equivalent operator.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Base.:==-Union{Tuple{T}, Tuple{T,T}} where T<:Hamiltonian.Prerequisites.Trees.AbstractTree",
    "page": "Trees",
    "title": "Base.:==",
    "category": "method",
    "text": "==(t1::T,t2::T) where T<:AbstractTree -> Bool\nisequal(t1::T,t2::T) where T<:AbstractTree -> Bool\n\nOverloaded equivalent operator.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Base.append!-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},AbstractTree{N,D}}} where D where N",
    "page": "Trees",
    "title": "Base.append!",
    "category": "method",
    "text": "append!(tree::AbstractTree{N,D},subtree::AbstractTree{N,D}) where {N,D} -> typeof(tree)\nappend!(tree::AbstractTree{N,D},node::Union{N,Nothing},subtree::AbstractTree{N,D}) where {N,D} -> typeof(tree)\n\nAppend a subtree to a tree.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Base.delete!-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},N}} where D where N",
    "page": "Trees",
    "title": "Base.delete!",
    "category": "method",
    "text": "delete!(tree::AbstractTree{N,D},node::N) where {N,D} -> typeof(tree)\n\nDelete a node and all its descendants from a tree.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Base.eltype-Tuple{Hamiltonian.Prerequisites.Trees.AbstractTree}",
    "page": "Trees",
    "title": "Base.eltype",
    "category": "method",
    "text": "eltype(tree::AbstractTree)\neltype(::Type{<:AbstractTree{N,D}}) where {N,D}\n\nGet the eltype of a tree.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Base.empty!-Tuple{Hamiltonian.Prerequisites.Trees.AbstractTree}",
    "page": "Trees",
    "title": "Base.empty!",
    "category": "method",
    "text": "empty!(tree::AbstractTree) -> typeof(tree)\n\nEmpty a tree.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Base.empty-Union{Tuple{AbstractTree{N,D}}, Tuple{D}, Tuple{N}} where D where N",
    "page": "Trees",
    "title": "Base.empty",
    "category": "method",
    "text": "empty(tree::AbstractTree)\n\nConstruct an empty tree of the same type with the input one.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Base.getindex-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},N}} where D where N",
    "page": "Trees",
    "title": "Base.getindex",
    "category": "method",
    "text": "getindex(tree::AbstractTree{N,D},node::N) where {N,D} -> N\n\nGet the data of a tree\'s node.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Base.haskey-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},N}} where D where N",
    "page": "Trees",
    "title": "Base.haskey",
    "category": "method",
    "text": "haskey(tree::AbstractTree{N,D},node::N) where {N,D} -> Bool\n\nCheck whether a node is in a tree.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Base.keys-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},TreeIteration}, Tuple{AbstractTree{N,D},TreeIteration,Union{Nothing, N}}} where D where N",
    "page": "Trees",
    "title": "Base.keys",
    "category": "method",
    "text": "keys(tree::AbstractTree{N,D},::TreeDepth,node::Union{N,Nothing}=tree|>root) where {N,D}\nkeys(tree::AbstractTree{N,D},::TreeWidth,node::Union{N,Nothing}=tree|>root) where {N,D}\n\nIterate over a tree\'s nodes starting from a certain node by depth first search or width first search.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Base.keytype-Tuple{Hamiltonian.Prerequisites.Trees.AbstractTree}",
    "page": "Trees",
    "title": "Base.keytype",
    "category": "method",
    "text": "keytype(tree::AbstractTree)\nkeytype(::Type{<:AbstractTree{N,D}}) where {N,D}\n\nGet a tree\'s node type.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Base.length-Tuple{Hamiltonian.Prerequisites.Trees.AbstractTree}",
    "page": "Trees",
    "title": "Base.length",
    "category": "method",
    "text": "length(tree::AbstractTree) -> Int\n\nGet the number of a tree\'s nodes.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Base.pairs-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree,TreeIteration}, Tuple{AbstractTree,TreeIteration,Union{Nothing, N}}} where D where N",
    "page": "Trees",
    "title": "Base.pairs",
    "category": "method",
    "text": "pairs(tree::AbstractTree,::TreeDepth,node::Union{N,Nothing}=tree|>root) where {N,D}\npairs(tree::AbstractTree,::TreeWidth,node::Union{N,Nothing}=tree|>root) where {N,D}\n\nIterate over a tree\'s (node,data) pairs starting from a certain node by depth first search or width first search.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Base.push!-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},N,D}} where D where N",
    "page": "Trees",
    "title": "Base.push!",
    "category": "method",
    "text": "push!(tree::AbstractTree{N,D},node::N,data::D) where {N,D} -> typeof(tree)\npush!(tree::AbstractTree{N,D},parent::Union{N,Nothing},node::N,data::D) where {N,D} -> typeof(tree)\n\nPush a new node to a tree. When parent is nothing, this function set the root node of an empty tree.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Base.setindex!-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree{N,D},D,N}} where D where N",
    "page": "Trees",
    "title": "Base.setindex!",
    "category": "method",
    "text": "setindex!(tree::AbstractTree{N,D},data::D,node::N) where {N,D}\n\nSet the data of a tree\'s node.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Base.valtype-Tuple{Hamiltonian.Prerequisites.Trees.AbstractTree}",
    "page": "Trees",
    "title": "Base.valtype",
    "category": "method",
    "text": "valtype(tree::AbstractTree)\nvaltype(::Type{<:AbstractTree{N,D}}) where {N,D}\n\nGet a tree\'s data type.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Base.values-Union{Tuple{D}, Tuple{N}, Tuple{AbstractTree,TreeIteration}, Tuple{AbstractTree,TreeIteration,Union{Nothing, N}}} where D where N",
    "page": "Trees",
    "title": "Base.values",
    "category": "method",
    "text": "values(tree::AbstractTree,::TreeDepth,node::Union{N,Nothing}=tree|>root) where {N,D}\nvalues(tree::AbstractTree,::TreeWidth,node::Union{N,Nothing}=tree|>root) where {N,D}\n\nIterate over a tree\'s data starting from a certain node by depth first search or width first search.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/Trees.html#Manual-1",
    "page": "Trees",
    "title": "Manual",
    "category": "section",
    "text": "Modules=[Trees]\nOrder=  [:module,:constant,:type,:macro,:function]"
},

{
    "location": "man/Prerequisites/NamedVectors.html#",
    "page": "Named vectors",
    "title": "Named vectors",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Prerequisites.NamedVectorspush!(LOAD_PATH,\"../../../../src/\")\nusing Hamiltonian.Prerequisites.NamedVectors"
},

{
    "location": "man/Prerequisites/NamedVectors.html#Named-vectors-1",
    "page": "Named vectors",
    "title": "Named vectors",
    "category": "section",
    "text": "A named vector is similiar to a named tuple, which associate each of its values with a name. Although the names of a named vector cannot be changed, the values can be modified if needed. In contrast to the predefined NamedTuple in Julia, which employs the names as type parameters, we just implement a named vector as a composite struct equipped with the getindex and setindex! functions, with the fieldnames being its names. This simple implementation makes it possible to define your own concrete named vector with any of your preferred type names, and ensures that all instances of a certain concrete named vector share the same names. Therefore, if you are familiar with Python, you will find that our named vector is more qualified to be the counterpart of the namedtuple in Python than the default Julia implementation. Last but not least important, it is also worth noted that a named vector is not a vector, as is similar to that a named tuple is not a tuple in Julia. This results from our basic expectation that a named vector should be more like a tuple other than a vector so that not all operations valid to vectors are also valid to named vectors."
},

{
    "location": "man/Prerequisites/NamedVectors.html#NamedVector-1",
    "page": "Named vectors",
    "title": "NamedVector",
    "category": "section",
    "text": "NamedVector defines the abstract type for all concrete named vectors.Main features include:Values can be accessed or modified either by the . operator or by the [] operator.\nComparisons, such as ≡, ≢, ==, ≠, >, <, ≥, ≤ are supported. Therefore a vector of named vectors can be sorted by the default sort function.\nHash is supported by hash. Therefore, a named vector can be used as the key of a dict or set.\nIteration over its fieldnames is supported by keys, over its values is supported by values, over its field-value pairs is supported by pairs. A reverse iteration is also supported.To subtype it, please note:A concrete type can be either mutable or immutable as you need, which is different from tuples.\nThe fields of a concrete type can be of the same type or not. For the former, we denote the named vector as \"homogeneous\" while for the latter as \"inhomogeneous\". For homogeneous ones, we define a sub abstract type, HomoNamedVector for further optimization of the default methods. See HomoNamedVector below.\nIt is recommended to overload the Base.fieldnames function for concrete subtypes to ensure type stability and improve efficiency, which though is not a necessity. A template for such an overloading is\nBase.fieldnames(Type{<:YourNamedVector})=(:fieldname1,:fieldname2,...)\nFor all concrete subtypes, if inner constructors are defined, the one which has the same interface with the default one must be implemented. Otherwise, some functionalities will not work.\nArithmetic operations, such as +, -, *, /, %, ÷, etc. are NOT supported. However, the function map is implemented, which can help users do the overloadings of these operations.We define a macro @namedvector as the type factory to decorate a \"raw\" struct to be a subtype of NamedVector. Here, \"raw\" means the struct to be decorated has no explicit supertype other than Any, neither inner constructors as well. For example,@namedvector mutable struct InHomoNV\n    scope::String\n    site::Int\nendThis macro encapsulate the overloading of Base.fieldnames, and you have no need to do this by hand any more."
},

{
    "location": "man/Prerequisites/NamedVectors.html#HomoNamedVector-1",
    "page": "Named vectors",
    "title": "HomoNamedVector",
    "category": "section",
    "text": "HomoNamedVector is the subtype of NamedVector that of all its fields share the same type. Compared to NamedVector, one more default method is implemented with HomoNamedVector, i.e. eltype, which returns the type of its fields. This function ensures the type stability of all the methods that involves an iteration of the field values of a named vector. Therefore, homogeneous named vector are usually more efficient than inhomogeneous ones. Use homogeneous ones as much as possible unless the code efficiency does not matter.To subtype HomoNamedVector, all the suggestions mentioned in the previous subsection for NamedVector also applies. A recommended template for a subtype is[mutable] struct YourNamedVector{T} <: HomoNamedVector{T}\n    filed1::T\n    filed2::T\n    ...\nendWe also provide a macro @homonamedvector to help the definition of concrete homogeneous named vector, where you only need specify the type name, field names, data type and optionally whether the subtype is mutable. For example,@homonamedvector HomoNVWithoutParameter (:scope,:site) Int mutable=true\n@homonamedvector HomoNVWithParameter (:scope,:site) (<:Real) mutable=trueThis macro also integrates the Base.fieldnames function, thus its overloading by hand is on longer needed."
},

{
    "location": "man/Prerequisites/NamedVectors.html#Hamiltonian.Prerequisites.NamedVectors.HomoNamedVector",
    "page": "Named vectors",
    "title": "Hamiltonian.Prerequisites.NamedVectors.HomoNamedVector",
    "category": "type",
    "text": "HomoNamedVector{T}\n\nAbstract type for all homogeneous named vectors.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/NamedVectors.html#Hamiltonian.Prerequisites.NamedVectors.NamedVector",
    "page": "Named vectors",
    "title": "Hamiltonian.Prerequisites.NamedVectors.NamedVector",
    "category": "type",
    "text": "NamedVector\n\nAbstract type for all named vectors.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/NamedVectors.html#Hamiltonian.Prerequisites.NamedVectors.@homonamedvector",
    "page": "Named vectors",
    "title": "Hamiltonian.Prerequisites.NamedVectors.@homonamedvector",
    "category": "macro",
    "text": "@homonamedvector typename fieldnames dtype::Union{Expr,Symbol}=:nothing mutable::Union{Expr,Bool}=false\n\nConstruct a concrete homogeneous named vector with the type name being typename and the fieldnames specified by fieldnames, and optionally, the type parameters specified by dtype.mutable can be used as a keyword argument to determine whether the concrete type is mutable.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/NamedVectors.html#Hamiltonian.Prerequisites.NamedVectors.@namedvector-Tuple{Expr}",
    "page": "Named vectors",
    "title": "Hamiltonian.Prerequisites.NamedVectors.@namedvector",
    "category": "macro",
    "text": "@namedvector structdef::Expr\n\nDecorate a \"raw\" struct to be a subtype of NamedVector. Here, \"raw\" means that the input struct has no explicit supertype and no inner constructors.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/NamedVectors.html#Base.:<-Tuple{Hamiltonian.Prerequisites.NamedVectors.NamedVector,Hamiltonian.Prerequisites.NamedVectors.NamedVector}",
    "page": "Named vectors",
    "title": "Base.:<",
    "category": "method",
    "text": "<(nv1::NamedVector,nv2::NamedVector) -> Bool\nisless(nv1::NamedVector,nv2::NamedVector) -> Bool\n\nCompare two named vectors and judge whether the first is less than the second.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/NamedVectors.html#Base.:==-Tuple{Hamiltonian.Prerequisites.NamedVectors.NamedVector,Hamiltonian.Prerequisites.NamedVectors.NamedVector}",
    "page": "Named vectors",
    "title": "Base.:==",
    "category": "method",
    "text": "==(nv1::NamedVector,nv2::NamedVector) -> Bool\nisequal(nv1::NamedVector,nv2::NamedVector) -> Bool\n\nOverloaded equivalent operator. Two named vector are equal to each other if and only if their keys as well as their values are equal to each other.\n\nnote: Note\nIt is not necessary for two named vectors to be of the same concrete type to be equal to each other.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/NamedVectors.html#Base.convert-Tuple{Type{Tuple},Hamiltonian.Prerequisites.NamedVectors.NamedVector}",
    "page": "Named vectors",
    "title": "Base.convert",
    "category": "method",
    "text": "convert(::Type{Tuple},nv::NamedVector) -> Tuple\nconvert(::Type{NV},nv::Tuple) where NV<:NamedVector -> NV\n\nConvert a named vector to tuple and vice versa.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/NamedVectors.html#Base.eltype-Union{Tuple{Type{#s68} where #s68<:HomoNamedVector{T}}, Tuple{T}} where T",
    "page": "Named vectors",
    "title": "Base.eltype",
    "category": "method",
    "text": "eltype(::Type{NV}) where NV<:HomoNamedVector{T} where T\neltype(nv::HomoNamedVector)\n\nGet the type parameter of a concrete HomoNamedVector.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/NamedVectors.html#Base.getindex-Tuple{Hamiltonian.Prerequisites.NamedVectors.NamedVector,Int64}",
    "page": "Named vectors",
    "title": "Base.getindex",
    "category": "method",
    "text": "getindex(nv::NamedVector,index::Int)\n\nGet the value by the [] syntax.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/NamedVectors.html#Base.hash-Tuple{Hamiltonian.Prerequisites.NamedVectors.NamedVector,UInt64}",
    "page": "Named vectors",
    "title": "Base.hash",
    "category": "method",
    "text": "hash(nv::NamedVector,h::UInt)\n\nHash a concrete NamedVector.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/NamedVectors.html#Base.iterate",
    "page": "Named vectors",
    "title": "Base.iterate",
    "category": "function",
    "text": "iterate(nv::NamedVector,state=1)\niterate(rv::Iterators.Reverse{<:NamedVector},state=length(rv.itr))\n\nIterate or reversely iterate over the values of a concrete NamedVector.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/NamedVectors.html#Base.keys-Tuple{Hamiltonian.Prerequisites.NamedVectors.NamedVector}",
    "page": "Named vectors",
    "title": "Base.keys",
    "category": "method",
    "text": "keys(nv::NamedVector) -> NTuple(nv|>fieldcount,Symbol)\nvalues(nv::NamedVector) -> Tuple\npairs(nv::NamedVector) -> Base.Generator\n\nIterate over the names.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/NamedVectors.html#Base.length-Union{Tuple{Type{NV}}, Tuple{NV}} where NV<:Hamiltonian.Prerequisites.NamedVectors.NamedVector",
    "page": "Named vectors",
    "title": "Base.length",
    "category": "method",
    "text": "length(::Type{NV}) where NV<:NamedVector -> Int\nlength(nv::NamedVector) -> Int\n\nGet the length of a concrete NamedVector.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/NamedVectors.html#Base.map-Union{Tuple{NV}, Tuple{Any,Vararg{NV,N} where N}} where NV<:Hamiltonian.Prerequisites.NamedVectors.NamedVector",
    "page": "Named vectors",
    "title": "Base.map",
    "category": "method",
    "text": "map(f,nvs::NV...) where NV<:NamedVector -> NV\n\nApply function f elementwise on the input named vectors.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/NamedVectors.html#Base.replace-Tuple{Hamiltonian.Prerequisites.NamedVectors.NamedVector}",
    "page": "Named vectors",
    "title": "Base.replace",
    "category": "method",
    "text": "replace(nv::NamedVector;kwargs...) -> typeof(nv)\n\nReturn a copy of a concrete NamedVector with some of the field values replaced by the keyword arguments.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/NamedVectors.html#Base.setindex!-Tuple{Hamiltonian.Prerequisites.NamedVectors.NamedVector,Any,Int64}",
    "page": "Named vectors",
    "title": "Base.setindex!",
    "category": "method",
    "text": "setindex!(nv::NamedVector,value,index::Int)\n\nSet the value by the [] syntax if mutable.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/NamedVectors.html#Base.show-Tuple{IO,Hamiltonian.Prerequisites.NamedVectors.NamedVector}",
    "page": "Named vectors",
    "title": "Base.show",
    "category": "method",
    "text": "show(io::IO,nv::NamedVector)\n\nShow a concrete NamedVector.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/NamedVectors.html#Base.zero-Union{Tuple{Type{NV}}, Tuple{NV}} where NV<:Hamiltonian.Prerequisites.NamedVectors.NamedVector",
    "page": "Named vectors",
    "title": "Base.zero",
    "category": "method",
    "text": "zero(::Type{NV}) where NV<:NamedVector -> NV\nzero(nv::NamedVector) -> typeof(nv)\n\nGet a concrete NamedVector with all values being zero.\n\n\n\n\n\n"
},

{
    "location": "man/Prerequisites/NamedVectors.html#Manual-1",
    "page": "Named vectors",
    "title": "Manual",
    "category": "section",
    "text": "Modules=[NamedVectors]\nOrder=  [:module,:constant,:type,:macro,:function]"
},

{
    "location": "man/Mathematics/Introduction.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Mathematics"
},

{
    "location": "man/Mathematics/Introduction.html#Introduction-1",
    "page": "Introduction",
    "title": "Introduction",
    "category": "section",
    "text": "This module contains the mathematical prerequisites of the Hamiltonian package.Pages=[\n    \"Combinatorics.md\",\n    \"VectorSpaces.md\",\n    \"AlgebraOverFields.md\",\n    \"QuantumNumbers.md\",\n    ]\nDepth=2"
},

{
    "location": "man/Mathematics/Combinatorics.html#",
    "page": "Combinatorics",
    "title": "Combinatorics",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Mathematics.Combinatorics"
},

{
    "location": "man/Mathematics/Combinatorics.html#Combinatorics-1",
    "page": "Combinatorics",
    "title": "Combinatorics",
    "category": "section",
    "text": "This module implements the combinations and permutations of an indexable object, with duplicate elements allowed or not. Compared to another Julia package Combinatorics, the iterators return tuples instead of vectors, which greatly decreases the momory allocation times and improves the code efficiency."
},

{
    "location": "man/Mathematics/Combinatorics.html#AbstractCombinatorics-1",
    "page": "Combinatorics",
    "title": "AbstractCombinatorics",
    "category": "section",
    "text": "AbstractCombinatorics{M,C} is the abstract type of all combinatoric algorithms. It has two type parameters:M: the number of elements to be taken\nC: the type of the collection of candidate elementsTo avoid momery allocation, the iteration of a concrete combinatoric algorithm returns a tuple, whose length is M and eltype is eltype(C)."
},

{
    "location": "man/Mathematics/Combinatorics.html#Combinations-and-DulCombinations-1",
    "page": "Combinatorics",
    "title": "Combinations and DulCombinations",
    "category": "section",
    "text": "Combinations{M,C} and DulCombinations generate all the combinations of M elements from an indexable collection whose type is C, with the differences being that the former forbids duplicate elements in the combinations while the latter allows."
},

{
    "location": "man/Mathematics/Combinatorics.html#Permutations-and-DulPermutations-1",
    "page": "Combinatorics",
    "title": "Permutations and DulPermutations",
    "category": "section",
    "text": "Permutations{M,C} and DulPermutations generate all the permutations of M elements from an indexable collection whose type is C, with the differences being that the former forbids duplicate elements in the permutations while the latter allows."
},

{
    "location": "man/Mathematics/Combinatorics.html#Hamiltonian.Mathematics.Combinatorics.AbstractCombinatorics",
    "page": "Combinatorics",
    "title": "Hamiltonian.Mathematics.Combinatorics.AbstractCombinatorics",
    "category": "type",
    "text": "AbstractCombinatorics{M,C}\n\nAbstract combinatoric algorithms.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/Combinatorics.html#Hamiltonian.Mathematics.Combinatorics.Combinations",
    "page": "Combinatorics",
    "title": "Hamiltonian.Mathematics.Combinatorics.Combinations",
    "category": "type",
    "text": "Combinations{M}(contents::C) where {M,C}\n\nCombinations of M elements from contents. Duplicates are not allowed.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/Combinatorics.html#Hamiltonian.Mathematics.Combinatorics.DulCombinations",
    "page": "Combinatorics",
    "title": "Hamiltonian.Mathematics.Combinatorics.DulCombinations",
    "category": "type",
    "text": "DulCombinations{M}(contents::C) where {M,C}\n\nCombinations of M elements from contents. Duplicates are allowed.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/Combinatorics.html#Hamiltonian.Mathematics.Combinatorics.DulPermutations",
    "page": "Combinatorics",
    "title": "Hamiltonian.Mathematics.Combinatorics.DulPermutations",
    "category": "type",
    "text": "DulPermutations{M}(contents::C) where {M,C}\n\nPermutations of M elements from contents. Duplicates are not allowed.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/Combinatorics.html#Hamiltonian.Mathematics.Combinatorics.Permutations",
    "page": "Combinatorics",
    "title": "Hamiltonian.Mathematics.Combinatorics.Permutations",
    "category": "type",
    "text": "Permutations{M}(contents::C) where {M,C}\n\nPermutations of M elements from contents. Duplicates are allowed.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/Combinatorics.html#Manul-1",
    "page": "Combinatorics",
    "title": "Manul",
    "category": "section",
    "text": "Modules=[Combinatorics]\nOrder=  [:module,:constant,:type,:macro,:function]"
},

{
    "location": "man/Mathematics/VectorSpaces.html#",
    "page": "Vector spaces",
    "title": "Vector spaces",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Mathematics.VectorSpaces"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Vector-spaces-1",
    "page": "Vector spaces",
    "title": "Vector spaces",
    "category": "section",
    "text": "A vector space is a linear space, in which the addition of vectors and multiplication of a vector by a scalar are defined.Vector spaces are frequently encountered in physics, e.g. the Hilbert space in quantum mechanics. In this submodule, we only implement those with finite dimensions. We want to remark that in our implementation, a vector space is a subtype of an abstract vector, therefore, the bases always possess a order, which means, two vector spaces are not considered to be equal to each other even if their corresponding actual mathmatical spaces are the same but the the orders of the bases are different."
},

{
    "location": "man/Mathematics/VectorSpaces.html#VectorSpace-1",
    "page": "Vector spaces",
    "title": "VectorSpace",
    "category": "section",
    "text": "VectorSpace{B} is the abstaction of a vector space, which has only one type parameter:B<:Any: the type of the bases of the vector spaceBasically, a subtype should implement the following 3 methods:dimension(vs::VectorSpace) -> Int\nGet the dimension of a vector space\nBase.getindex(vs::VectorSpace{B},i::Int)  where B -> B\nGet the ith basis of a vector space\nBase.searchsortedfirst(vs::VectorSpace{B},basis::B) where B -> Int\nSearch the index of a basis in a vector spaceHowever, we provide several interfaces, including type traits and methods to deal with common situations:A vector space whose bases are stored in a table under the attribute name :table can be ascribed to the HasTable trait and the TableSorted trait. Specifically, the first trait must be implemented as\nHasTable(::Type{SubType})=HasTable(true)\nWhile, if the table is unsorted, the second trait should be implemented as\nTableSorted(::Type{SubType})=TableSorted(false)\nand if the table is sorted, the second trait should be implemented as\nTableSorted(::Type{SubType})=TableSorted(true)\nA vector space whose bases may be represented by a multiindex (Cartesian index) can be ascribed to the traits IsMultiIndexable and MultiIndexOrderStyle. Specifically, the first trait must be implemented as\nIsMultiIndexable(::Type{SubType})=IsMultiIndexable(true)\nWhile, if the order style of the multiindex is C/C++ like, the second trait shoule be implemented as\nMultiIndexOrderStyle(::Type{SubType})=MultiIndexOrderStyle(\'C\')\nand if the order style is Fortran like, the second trait shoule be implemented as\nMultiIndexOrderStyle(::Type{SubType})=MultiIndexOrderStyle(\'F\')\nFurthermore, it should implement the following methods\nrank(::Type{SubType}) -> Int\nGet the rank of a multiindexable vector space.\ndims(vs::SubType) -> NTuple{vs|>typeof|>rank,Int}\nGet the dimensions along each axes of a multiindexable vector space.\ninds(basis,vs::SubType) ->  NTuple{vs|>typeof|>rank,Int}\nGet the Cartesian index representation of a basis in a multiindexable vector space.\neltype(SubType).name.wrapper(index::NTuple{N,Int},vs::SubType)\nConstruct a basis from a tuple of integers and a multiindexable vector space.\nNote that a multiindexable vector space can also have a sorted or unsorted table. But then the trait MultiIndexOrderStyle takes no effects and the sequences of its bases will be completely determined by its attribute :table.If the type taits and methods are defined properly as stated above, the dimension, getindex and searchsortedfirst methods get default implementations. No need to worry about them any more.Other features includecomparison: == and isequal\niteration: iterate\ninquiry: size, findfirst and in"
},

{
    "location": "man/Mathematics/VectorSpaces.html#DirectVectorSpace-1",
    "page": "Vector spaces",
    "title": "DirectVectorSpace",
    "category": "section",
    "text": "DirectVectorSpace{S,B,N} is the simplest vector space, whose bases are stored in the attribute :table as an ntuple.The :table attribute can be sorted or unsorted, which is determined by the type parameter S, with \'T\' for sorted and \'F\' for unsorted."
},

{
    "location": "man/Mathematics/VectorSpaces.html#OrderedIndices-1",
    "page": "Vector spaces",
    "title": "OrderedIndices",
    "category": "section",
    "text": "OrderedIndices{N} defines the simplest abstract class of multiindexable vector spaces, whose bases are just tuples of integers.This class of vector spaces must have the following attribute: dims::NTuple{N,Int}: the dimesnions of the Cartesian indices along all axes"
},

{
    "location": "man/Mathematics/VectorSpaces.html#DirectIndices-1",
    "page": "Vector spaces",
    "title": "DirectIndices",
    "category": "section",
    "text": "DirectIndices{M,N} is the direct ordered Cartesian indices.It is worth noting thatIt can be C/C++ ordered or Fortran ordered depending on the first type parameter M, with \'C\' for the former and \'F\' the latter.\nFor its bases (Cartesian indices), there is no restriction except that they should be in the proper range defined by its dims."
},

{
    "location": "man/Mathematics/VectorSpaces.html#TabledIndices-1",
    "page": "Vector spaces",
    "title": "TabledIndices",
    "category": "section",
    "text": "TabledIndices{S,N} defines the tabled ordered Cartesian indices.Compared to DirectIndices, the bases of this kind of vector spaces are stored in the attribute :table, which must be a vector of tuple of integers. The :table attribute can be sorted or unsorted, which is determined by the type parameter S, with \'T\' for sorted and \'F\' for unsorted. This type suits the situations when the Cartesian indices are restricted by extra conditions except that propoesed by the attribute :dims."
},

{
    "location": "man/Mathematics/VectorSpaces.html#GradedVectorSpace-1",
    "page": "Vector spaces",
    "title": "GradedVectorSpace",
    "category": "section",
    "text": "GradedVectorSpace{G,B,V,T} defines the abstract type of graded vector spaces, which are vector spaces that have the extra structure of a grading, which is a decomposition of the vector space into a direct sum of vector subspaces.It has 4 type parametersG: the type of the grades\nB: the eltype of the subspaces\nV<:VectorSpace: the type of the subspaces\nT<:GradedTables{G,V}: the type of the subspaces\' contentsConcrete subtypes must have the following attribute::tables::T: the contents of the subspaces, which must be a GradedTables.Specifically, the dimension, getindex and searchsortedfirst methods are overloaded in support of various purposes. For details, please refer to the manual."
},

{
    "location": "man/Mathematics/VectorSpaces.html#Hamiltonian.Mathematics.VectorSpaces.DirectIndices",
    "page": "Vector spaces",
    "title": "Hamiltonian.Mathematics.VectorSpaces.DirectIndices",
    "category": "type",
    "text": "DirectIndices{M}(dims::NTuple{N,Int}) where {M,N}\n\nDirect ordered Cartesian indices.\n\nnote: Note\nIt can be C/C++ ordered or Fortran ordered depending on the first type parameter M, with \'C\' for the former and \'F\' the latter.\nFor its bases (Cartesian indices), there is no restriction except that they should be in the proper range defined by its dims.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Hamiltonian.Mathematics.VectorSpaces.DirectVectorSpace",
    "page": "Vector spaces",
    "title": "Hamiltonian.Mathematics.VectorSpaces.DirectVectorSpace",
    "category": "type",
    "text": "DirectVectorSpace{S}(table::NTuple{N,B}) where {S,B,N}\nDirectVectorSpace{S}(table...) where S\n\nSimplest vector space, whose bases are stored in the attribute :table as an ntuple.\n\nThe :table attribute can be sorted or unsorted, which is determined by the type parameter S, with \'T\' for sorted and \'F\' for unsorted.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Hamiltonian.Mathematics.VectorSpaces.GradedTables",
    "page": "Vector spaces",
    "title": "Hamiltonian.Mathematics.VectorSpaces.GradedTables",
    "category": "type",
    "text": "GradedTables(vs::Tuple,ks::Tuple)\nGradedTables(::Type{M},n::Int,gs::Val{GS}) where {M<:AbstractCombinatorics,GS}\n\nThe tables of a graded vector space.\n\nAlias for Base.Iterators.Pairs{G,V,KS<:Tuple,VS<:Tuple}.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Hamiltonian.Mathematics.VectorSpaces.GradedVectorSpace",
    "page": "Vector spaces",
    "title": "Hamiltonian.Mathematics.VectorSpaces.GradedVectorSpace",
    "category": "type",
    "text": "GradedVectorSpace{G,B,V<:VectorSpace,T<:GradedTables{G,V}} <: VectorSpace{Tuple{G,B}}\n\nAbstract type of graded vector spaces.\n\nA graded vector space is a vector space that has the extra structure of a grading, which is a decomposition of the vector space into a direct sum of vector subspaces.\n\nConcrete subtypes must have the following attribute:\n\n:tables::GradedTables{G,V} where {G,V<:VectorSpace}: the tables of the subspaces.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Hamiltonian.Mathematics.VectorSpaces.HasTable",
    "page": "Vector spaces",
    "title": "Hamiltonian.Mathematics.VectorSpaces.HasTable",
    "category": "type",
    "text": "HasTable(B::Bool)\nHasTable(::Type{<:VectorSpace})\n\nTrait of whether a subtype of VectorSpace has the attribute :table.\n\nOnly two instances are allowed, the first of which is the default for a subtype:\n\nHasTable(false): indication of not having the attribute :table\nHasTable(true): indication of having the attribute :table\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Hamiltonian.Mathematics.VectorSpaces.IsMultiIndexable",
    "page": "Vector spaces",
    "title": "Hamiltonian.Mathematics.VectorSpaces.IsMultiIndexable",
    "category": "type",
    "text": "IsMultiIndexable(B::Bool)\nIsMultiIndexable(::Type{<:VectorSpace})\n\nTrait of whether the bases of a subtype of VectorSpace can be represented by multiindices (Cartesian indices).\n\nOnly two instances are allowed, the first of which is the default for a subtype:\n\nIsMultiIndexable(false): indication of irrepresentability by Cartesian indices\nIsMultiIndexable(true): indication of representability by Cartesian indices\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Hamiltonian.Mathematics.VectorSpaces.MultiIndexOrderStyle",
    "page": "Vector spaces",
    "title": "Hamiltonian.Mathematics.VectorSpaces.MultiIndexOrderStyle",
    "category": "type",
    "text": "MultiIndexOrderStyle(M::Char)\nMultiIndexOrderStyle(::Type{<:VectorSpace})\n\nTrait of the order style of the Cartesian-index representation of the bases of a subtype of VectorSpace.\n\nOnly two instances are allowed, the first of which is the default for a subtype:\n\nMultiIndexOrderStyle(C): indication of C/C++ order style\nMultiIndexOrderStyle(F): indication of Fortran order style\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Hamiltonian.Mathematics.VectorSpaces.OrderedIndices",
    "page": "Vector spaces",
    "title": "Hamiltonian.Mathematics.VectorSpaces.OrderedIndices",
    "category": "type",
    "text": "OrderedIndices{N} <: VectorSpace{NTuple{N,Int}}\n\nThe simplest abstract class of multiindexable vector spaces, whose bases are just tuples of integers.\n\nThis class of vector spaces must have the following attribute: dims::NTuple{N,Int}: the dimesnions of the Cartesian indices along all axes\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Hamiltonian.Mathematics.VectorSpaces.TableSorted",
    "page": "Vector spaces",
    "title": "Hamiltonian.Mathematics.VectorSpaces.TableSorted",
    "category": "type",
    "text": "TableSorted(B::Bool)\nTableSorted(::Type{<:VectorSpace})\n\nTrait of whether the attribute :table of a subtype of VectorSpace is sorted.\n\nOnly two instances are allowed, the first of which is the default for a subtype:\n\nTableSorted(false): indication of unsorted attribute :table\nTableSorted(true): indication of sorted attribute :table\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Hamiltonian.Mathematics.VectorSpaces.TabledIndices",
    "page": "Vector spaces",
    "title": "Hamiltonian.Mathematics.VectorSpaces.TabledIndices",
    "category": "type",
    "text": "TabledIndices{S}(dims::NTuple{N,Int},table::Vector{NTuple{N,Int}}) where {S,N}\nTabledIndices{N}(::Type{M},len::Int) where {N,M<:AbstractCombinatorics}\n\nTabled ordered Cartesian indices.\n\nCompared to DirectIndices, the bases of this kind of vector spaces are stored in the attribute :table, which must be a vector of tuple of integers. The :table attribute can be sorted or unsorted, which is determined by the type parameter S, with \'T\' for sorted and \'F\' for unsorted. This type suits the situations when the Cartesian indices are restricted by extra conditions except that propoesed by the attribute :dims.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Hamiltonian.Mathematics.VectorSpaces.VectorSpace",
    "page": "Vector spaces",
    "title": "Hamiltonian.Mathematics.VectorSpaces.VectorSpace",
    "category": "type",
    "text": "VectorSpace{B} <: AbstractVector{B}\n\nAbstract vector space.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Hamiltonian.Prerequisites.Interfaces.:⊕-Union{Tuple{B}, Tuple{B,B}} where B",
    "page": "Vector spaces",
    "title": "Hamiltonian.Prerequisites.Interfaces.:⊕",
    "category": "method",
    "text": "⊕(basis1::B,basis2::B) where B -> DirectVectorSpace{\'F\',B,2}\n⊕(basis::B,vs::DirectVectorSpace{<:Any,B}) where B -> DirectVectorSpace{\'F\',B}\n⊕(vs::DirectVectorSpace{<:Any,B},basis::B) where B -> DirectVectorSpace{\'F\',B}\n⊕(vs1::DirectVectorSpace{<:Any,B},vs2::DirectVectorSpace{<:Any,B}) where B -> DirectVectorSpace{\'F\',B}\n\nGet the direct sum between bases or direct vector spaces.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Hamiltonian.Prerequisites.Interfaces.degree-Union{Tuple{G}, Tuple{G,GradedVectorSpace{G,B,V,T} where T<:(Pairs{G,V,KS,VS} where VS<:Tuple where KS<:Tuple) where V<:VectorSpace where B}} where G",
    "page": "Vector spaces",
    "title": "Hamiltonian.Prerequisites.Interfaces.degree",
    "category": "method",
    "text": "degree(g::G,vs::GradedVectorSpace{G}) where G -> Int\n\nGet the degree of a vector subspace whose grade are represented by g.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Hamiltonian.Prerequisites.Interfaces.dimension-Tuple{Hamiltonian.Mathematics.VectorSpaces.GradedVectorSpace}",
    "page": "Vector spaces",
    "title": "Hamiltonian.Prerequisites.Interfaces.dimension",
    "category": "method",
    "text": "dimension(vs::GradedVectorSpace) -> Int\ndimension(vs::GradedVectorSpace{G},g::G) where G -> Int\ndimension(vs::GradedVectorSpace{G},gs::NTuple{N,G}) where {G,N} -> Int\n\nGet the dimension of the whole graded vector space or some vector subspaces.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Hamiltonian.Prerequisites.Interfaces.dimension-Tuple{Hamiltonian.Mathematics.VectorSpaces.VectorSpace}",
    "page": "Vector spaces",
    "title": "Hamiltonian.Prerequisites.Interfaces.dimension",
    "category": "method",
    "text": "dimension(vs::VectorSpace) -> Int\n\nThe dimension of a vector space.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Hamiltonian.Prerequisites.Interfaces.rank-Tuple{Base.Iterators.Pairs{G,V,KS,VS} where VS<:Tuple where KS<:Tuple where V where G}",
    "page": "Vector spaces",
    "title": "Hamiltonian.Prerequisites.Interfaces.rank",
    "category": "method",
    "text": "rank(tables::GradedTables) -> Int\nrank(::Type{T}) where {T<:GradedTables} -> Int\n\nGet the total number of keys or values of a graded tables.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Hamiltonian.Prerequisites.Interfaces.rank-Tuple{Hamiltonian.Mathematics.VectorSpaces.GradedVectorSpace}",
    "page": "Vector spaces",
    "title": "Hamiltonian.Prerequisites.Interfaces.rank",
    "category": "method",
    "text": "rank(vs::GradedVectorSpace) -> Int\nrank(::Type{V}) where {V<:GradedVectorSpace} -> Int\n\nGet the rank, i.e. the total number of vector subspaces.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Base.Sort.searchsortedfirst-Union{Tuple{B}, Tuple{G}, Tuple{GradedVectorSpace{G,B,V,T} where T<:(Pairs{G,V,KS,VS} where VS<:Tuple where KS<:Tuple) where V<:VectorSpace,Tuple{G,B}}} where B where G",
    "page": "Vector spaces",
    "title": "Base.Sort.searchsortedfirst",
    "category": "method",
    "text": "searchsortedfirst(vs::GradedVectorSpace{G,B},pair::Tuple{G,B}) where {G,B} -> Int\n\nFind the index of a grade-basis pair in a graded vector space.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Base.Sort.searchsortedfirst-Union{Tuple{B}, Tuple{VectorSpace{B},B}} where B",
    "page": "Vector spaces",
    "title": "Base.Sort.searchsortedfirst",
    "category": "method",
    "text": "searchsortedfirst(vs::VectorSpace{B},basis::B) where B -> Int\nsearchsortedfirst(vs,basis) -> Int\n\nSearch the index of a basis in a vector space.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Base.eltype-Tuple{Hamiltonian.Mathematics.VectorSpaces.GradedVectorSpace,Int64}",
    "page": "Vector spaces",
    "title": "Base.eltype",
    "category": "method",
    "text": "eltype(vs::GradedVectorSpace,g::Int)\neltype(::Type{V},g::Int) where {V<:GradedVectorSpace}\n\nGet the gth eltype of a graded vector space.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Base.findfirst-Union{Tuple{B}, Tuple{B,VectorSpace{B}}} where B",
    "page": "Vector spaces",
    "title": "Base.findfirst",
    "category": "method",
    "text": "findfirst(basis::B,vs::VectorSpace{B}) where B -> Int\nfindfirst(bases,vs::VectorSpace) -> NTuple{length(bases),Int}\n\nGet the index of a basis or the indices of a couple of bases in a vector space.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Base.getindex-Tuple{Hamiltonian.Mathematics.VectorSpaces.VectorSpace,Int64}",
    "page": "Vector spaces",
    "title": "Base.getindex",
    "category": "method",
    "text": "getindex(vs::VectorSpace,i::Int)\n\nGet the ith basis of a vector space by the [] operator.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Base.getindex-Union{Tuple{B}, Tuple{GradedVectorSpace{B,B1,V,T} where T<:(Pairs{B,V,KS,VS} where VS<:Tuple where KS<:Tuple) where V<:VectorSpace where B1,Tuple{B,Int64}}} where B",
    "page": "Vector spaces",
    "title": "Base.getindex",
    "category": "method",
    "text": "getindex(vs::GradedVectorSpace{B},pair::Tuple{B,Int}) where B\ngetindex(vs::GradedVectorSpace,i::Int)\n\nGet the basis of a graded vector space by a grade-index pair or by an index.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Base.in-Union{Tuple{B}, Tuple{B,VectorSpace{B}}} where B",
    "page": "Vector spaces",
    "title": "Base.in",
    "category": "method",
    "text": "in(basis::B,vs::VectorSpace{B}) where B -> Bool\n\nJudge whether a basis belongs to a vector space.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Base.iterate",
    "page": "Vector spaces",
    "title": "Base.iterate",
    "category": "function",
    "text": "iterate(vs::GradedVectorSpace,state=(1,1))\n\nIterate over the whole graded vector space.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Base.keys-Tuple{Hamiltonian.Mathematics.VectorSpaces.GradedVectorSpace}",
    "page": "Vector spaces",
    "title": "Base.keys",
    "category": "method",
    "text": "keys(vs::GradedVectorSpace) -> Tuple\nvalues(vs::GradedVectorSpace) -> Tuple\npairs(vs::GradedVectorSpace) -> GradedTables\n\nIterate over the keys, values or pairs of a graded vector space.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Base.keytype-Tuple{Base.Iterators.Pairs{G,V,KS,VS} where VS<:Tuple where KS<:Tuple where V where G,Int64}",
    "page": "Vector spaces",
    "title": "Base.keytype",
    "category": "method",
    "text": "keytype(tables::GradedTables,g::Int)\nkeytype(::Type{T},g::Int) where {T<:GradedTables}\n\nGet the gth keytype of a graded tables.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Base.keytype-Tuple{Hamiltonian.Mathematics.VectorSpaces.GradedVectorSpace,Int64}",
    "page": "Vector spaces",
    "title": "Base.keytype",
    "category": "method",
    "text": "keytype(vs::GradedVectorSpace,g::Int)\nkeytype(::Type{V},g::Int) where {V<:GradedVectorSpace}\n\nGet the gth keytype of a graded vector space.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Base.valtype-Tuple{Base.Iterators.Pairs{G,V,KS,VS} where VS<:Tuple where KS<:Tuple where V where G,Int64}",
    "page": "Vector spaces",
    "title": "Base.valtype",
    "category": "method",
    "text": "valtype(tables::GradedTables,g::Int)\nvaltype(::Type{T},g::Int) where {T<:GradedTables}\n\nGet the gth valtype of a graded tables.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Base.valtype-Tuple{Hamiltonian.Mathematics.VectorSpaces.GradedVectorSpace,Int64}",
    "page": "Vector spaces",
    "title": "Base.valtype",
    "category": "method",
    "text": "valtype(vs::GradedVectorSpace,g::Int)\nvaltype(::Type{V},g::Int) where {V<:GradedVectorSpace}\n\nGet the gth valtype of a graded vector space.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/VectorSpaces.html#Manul-1",
    "page": "Vector spaces",
    "title": "Manul",
    "category": "section",
    "text": "Modules=[VectorSpaces]\nOrder=  [:module,:constant,:type,:macro,:function]"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#",
    "page": "Algebra over fields",
    "title": "Algebra over fields",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Mathematics.AlgebraOverFieldspush!(LOAD_PATH,\"../../../../src/\")\nusing Hamiltonian.Mathematics.AlgebraOverFields"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Algebra-over-fields-1",
    "page": "Algebra over fields",
    "title": "Algebra over fields",
    "category": "section",
    "text": "An algebra over a field is a vector space over that field, in which a bilinear operator (often called the \"multiplication\") between vectors is defined.With the help of the structure constants of the algebra, the result of the bilinear operation between any arbitary two vectors can be expressed by a sum of individual ones. Therefore, in principle, an algebra can be represented by the complete basis set of its corresponding vector space and a rank-3 tensor encapsulating its structure constants. It is noted that the \"bilinear operation\" is not restriced to the usual multiplication only. For example, it is the commutator, which is a composition of the usual multiplication and subtraction (for any A and B, the commutator [A,B] is defined as [A,B]≝AB-BA) that serves as the bilinear operator for Lie algebras. In this module, for scalars in the field and elements in the algebra, we only provide the interfaces of the scalar multiplication (including the scalar division) bwteen a sclar and an element, the addition (including the subtraction) and the usual multiplication between two elements. Other complicated operations should be composed from these basic ones."
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#SimpleID-and-ID-1",
    "page": "Algebra over fields",
    "title": "SimpleID and ID",
    "category": "section",
    "text": "SimpleID is the building block of the id system of an algebra over a field, while ID defines the specific identifier of an element in that algebra.Generally, the usual multiplication between two elements of an algebra is not commutable, and the rank of the multiplication is just the add-up before the simplication with the help of the algebra structure. We thus use a simple id to denote a single basis of the corresponding vector space, and an id to denote the identifier of an element. With the help of the direct product (⊗) of two ids, an over complete id system designed for the whole algebra is constructed. This id system is redundant because it does not reflects the structure constants of the algebra, which reduces independent basis elements. Extra mechanisms should be provided to kill this redundancy, which goes beyond the current module. Users should define them themselves."
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#IdSpace-1",
    "page": "Algebra over fields",
    "title": "IdSpace",
    "category": "section",
    "text": "IdSpace defines the complete vector space that corresponds to an algebra.An id space uses a set of simple id bases to generated all ranked ids, with the former represented by an instance of DirectVectorSpace and stored in the attribute :sids, while the latter by an instance of GradedTables in the attribute :tables. For the sick of the orderings among different ids, an id will be represented by a multiindex, which are determined by the corresponding sequences of its simple ids in the basis set. The id system can be complete or over complete, depending on whether the attribute :tables contain multiindices that actually represent the same basis of the algebra. Besides, the methodsBase.getindex(idspace::IdSpace,i::Int) -> IDandBase.findfirst(id::ID,idspace::IdSpace) -> Intare provided to get the ith id and to find the index of an id of a idspace, respectively."
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Element-and-Elements-1",
    "page": "Algebra over fields",
    "title": "Element and Elements",
    "category": "section",
    "text": "Element defines a single element of an algebra while Elements defines an exprssion composed of several elements from an algebra.The first and second attributes of an Element must bevalue::Number: the coefficient of the element\nid::ID: the id of the elementArithmetic operations (+, -, *, /) bwteen a scalar, an Element or an Elements is defined. See Manual for details."
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Hamiltonian.Mathematics.AlgebraOverFields.Element",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Mathematics.AlgebraOverFields.Element",
    "category": "type",
    "text": "Element{N,V<:Number,I<:ID{<:NTuple{N,SimpleID}}}\n\nAn element of an algebra over a field.\n\nThe first and second attributes of an element must be\n\nvalue::Number: the coefficient of the element\nid::ID: the id of the element\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Hamiltonian.Mathematics.AlgebraOverFields.Elements",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Mathematics.AlgebraOverFields.Elements",
    "category": "type",
    "text": "Elements{I<:ID,M<:Element} <: AbstractDict{I,M}\n\nAn set of elements of an algebra over a field.\n\nAlias for Dict{I<:ID,M<:Element}. Similar iterms are automatically merged thanks to the id system.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Hamiltonian.Mathematics.AlgebraOverFields.Elements-Tuple{Any}",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Mathematics.AlgebraOverFields.Elements",
    "category": "method",
    "text": "Elements(ms)\nElements(ms::Pair{I,M}...) where {I<:ID,M<:Element}\nElements(ms::Element...)\n\nGet the set of elements with similar items merged.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Hamiltonian.Mathematics.AlgebraOverFields.ID",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Mathematics.AlgebraOverFields.ID",
    "category": "type",
    "text": "ID(ids::NTuple{N,SimpleID}) where N\nID(ids::SimpleID...)\nID(::Type{SID},attrs::Vararg{NTuple{N},M}) where {SID<:SimpleID,N,M}\n\nThe id system of an algebra over a field.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Hamiltonian.Mathematics.AlgebraOverFields.IdSpace",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Mathematics.AlgebraOverFields.IdSpace",
    "category": "type",
    "text": "IdSpace(sids::DirectVectorSpace,tables::GradedTables)\nIdSpace(::Type{M},sids::DirectVectorSpace,gs::Val{GS}) where {M<:AbstractCombinatorics,GS}\n\nThe graded id space for an algebra generated by a couple of basic simple ids.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Hamiltonian.Mathematics.AlgebraOverFields.SimpleID",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Mathematics.AlgebraOverFields.SimpleID",
    "category": "type",
    "text": "SimpleID <: NamedVector\n\nA simple id is the building block of the id system of an algebra over a field.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Hamiltonian.Mathematics.AlgebraOverFields.idtype-Union{Tuple{Type{#s15} where #s15<:(Element{N,#s14,I} where #s14<:Number)}, Tuple{I}, Tuple{N}} where I<:(Hamiltonian.Mathematics.AlgebraOverFields.ID{#s16} where #s16<:Tuple{Vararg{Hamiltonian.Mathematics.AlgebraOverFields.SimpleID,N}}) where N",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Mathematics.AlgebraOverFields.idtype",
    "category": "method",
    "text": "idtype(::Type{<:Element{N,<:Number,I}}) where {N,I<:ID{<:NTuple{N,SimpleID}}}\nidtype(m::Element)\n\nThe type of the id of an element.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Hamiltonian.Prerequisites.Interfaces.:⊗-Tuple{Hamiltonian.Mathematics.AlgebraOverFields.Element,Dict{I,M} where M<:Hamiltonian.Mathematics.AlgebraOverFields.Element where I<:Hamiltonian.Mathematics.AlgebraOverFields.ID}",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Prerequisites.Interfaces.:⊗",
    "category": "method",
    "text": "⊗(m::Element,ms::Elements) -> Elements\n⊗(ms::Elements,m::Element) -> Elements\n⊗(ms1::Elements,ms2::Elements) -> Elements\n\nOverloaded ⊗ operator for element-element multiplications of an algebra over a field.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Hamiltonian.Prerequisites.Interfaces.add!-Tuple{Dict{I,M} where M<:Hamiltonian.Mathematics.AlgebraOverFields.Element where I<:Hamiltonian.Mathematics.AlgebraOverFields.ID}",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Prerequisites.Interfaces.add!",
    "category": "method",
    "text": "add!(ms::Elements) -> typeof(ms)\nadd!(ms::Elements,::Tuple{}) -> typeof(ms)\nadd!(ms::Elements,m::Element) -> typeof(ms)\nadd!(ms::Elements,mms::Elements) -> typeof(ms)\n\nGet the inplace addition of elements to a set.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Hamiltonian.Prerequisites.Interfaces.rank-Union{Tuple{Type{#s68} where #s68<:(Element{N,V,I} where I<:(ID{#s69} where #s69<:Tuple{Vararg{SimpleID,N}}) where V<:Number)}, Tuple{N}} where N",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Prerequisites.Interfaces.rank",
    "category": "method",
    "text": "rank(::Type{<:Element{N}}) where N -> Int\nrank(m::Element) -> Int\n\nGet the rank of an element.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Hamiltonian.Prerequisites.Interfaces.rank-Union{Tuple{Type{#s68} where #s68<:ID{T}}, Tuple{T}} where T<:Tuple{Vararg{Hamiltonian.Mathematics.AlgebraOverFields.SimpleID,N} where N}",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Prerequisites.Interfaces.rank",
    "category": "method",
    "text": "rank(::Type{<:ID{T}}) where T<:Tuple{Vararg{SimpleID}} -> Int\nrank(id::ID) -> Int\n\nGet the rank of a composite id.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Hamiltonian.Prerequisites.Interfaces.sub!-Tuple{Dict{I,M} where M<:Hamiltonian.Mathematics.AlgebraOverFields.Element where I<:Hamiltonian.Mathematics.AlgebraOverFields.ID}",
    "page": "Algebra over fields",
    "title": "Hamiltonian.Prerequisites.Interfaces.sub!",
    "category": "method",
    "text": "sub!(ms::Elements) -> typeof(ms)\nsub!(ms::Elements,::Tuple{}) -> typeof(ms)\nsub!(ms::Elements,m::Element) -> typeof(ms)\nsub!(ms::Elements,mms::Elements) -> typeof(ms)\n\nGet the inplace subtraction of elements from a set.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Base.:*-Tuple{Hamiltonian.Mathematics.AlgebraOverFields.SimpleID,Hamiltonian.Mathematics.AlgebraOverFields.SimpleID}",
    "page": "Algebra over fields",
    "title": "Base.:*",
    "category": "method",
    "text": "*(sid1::SimpleID,sid2::SimpleID) -> ID\n*(sid::SimpleID,cid::ID) -> ID\n*(cid::ID,sid::SimpleID) -> ID\n*(cid1::ID,cid2::ID) -> ID\n\nGet the product of the id system.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Base.:*-Tuple{Number,Hamiltonian.Mathematics.AlgebraOverFields.Element}",
    "page": "Algebra over fields",
    "title": "Base.:*",
    "category": "method",
    "text": "*(factor::Number,m::Element) -> Element\n*(m::Element,factor::Number) -> Element\n*(m1::Element,m2::Element) -> Element\n*(factor::Number,ms::Elements) -> Elements\n*(ms::Elements,factor::Number) -> Elements\n*(m::Element,ms::Elements) -> Elements\n*(ms::Elements,m::Element) -> Elements\n*(ms1::Elements,ms2::Elements) -> Elements\n\nOverloaded * operator for element-scalar multiplications and element-element multiplications of an algebra over a field.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Base.:+-Tuple{Hamiltonian.Mathematics.AlgebraOverFields.Element}",
    "page": "Algebra over fields",
    "title": "Base.:+",
    "category": "method",
    "text": "+(m::Element) -> typeof(m)\n+(ms::Elements) -> typeof(ms)\n+(m::Element,::Tuple{}) -> typeof(m)\n+(::Tuple{},m::Element) -> typeof(m)\n+(ms::Elements,::Tuple{}) -> typeof(ms)\n+(::Tuple{},ms::Elements) -> typeof(ms)\n+(ms::Elements,m::Element) -> Elements\n+(m1::Element,m2::Element) -> Elements\n+(m::Element,ms::Elements) -> Elements\n+(ms1::Elements,ms2::Elements) -> Elements\n\nOverloaded + operator between elements of an algebra over a field.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Base.:--Tuple{Hamiltonian.Mathematics.AlgebraOverFields.Element}",
    "page": "Algebra over fields",
    "title": "Base.:-",
    "category": "method",
    "text": "-(m::Element) -> typeof(m)\n-(ms::Elements) -> typeof(ms)\n-(m::Element,::Tuple{}) -> typeof(m)\n-(::Tuple{},m::Element) -> typeof(m)\n-(ms::Elements,::Tuple{}) -> typeof(ms)\n-(::Tuple{},ms::Elements) -> typeof(ms)\n-(m1::Element,m2::Element) -> Elements\n-(m::Element,ms::Elements) -> Elements\n-(ms::Elements,m::Element) -> Elements\n-(ms1::Elements,ms2::Elements) -> Elements\n\nOverloaded - operator between elements of an algebra over a field.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Base.:/-Tuple{Hamiltonian.Mathematics.AlgebraOverFields.Element,Number}",
    "page": "Algebra over fields",
    "title": "Base.:/",
    "category": "method",
    "text": "/(m::Element,factor::Number)\n/(ms::Elements,factor::Number)\n\nOverloaded / operator for element-sclar division of an algebra over a field.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Base.:==-Tuple{Hamiltonian.Mathematics.AlgebraOverFields.Element,Hamiltonian.Mathematics.AlgebraOverFields.Element}",
    "page": "Algebra over fields",
    "title": "Base.:==",
    "category": "method",
    "text": "(m1::Element,m2::Element) -> Bool\nisequal(m1::Element,m2::Element) -> Bool\n\nCompare two elements and judge whether they are equal to each other.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Base.:^-Tuple{Hamiltonian.Mathematics.AlgebraOverFields.Element,Int64}",
    "page": "Algebra over fields",
    "title": "Base.:^",
    "category": "method",
    "text": "^(m::Element,n::Int)\n^(ms::Elements,n::Int)\n\nOverloaded ^ operator for element-integer power of an algebra over a field.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Base.empty-Tuple{Type{#s68} where #s68<:(Dict{I,M} where M<:Hamiltonian.Mathematics.AlgebraOverFields.Element where I<:Hamiltonian.Mathematics.AlgebraOverFields.ID)}",
    "page": "Algebra over fields",
    "title": "Base.empty",
    "category": "method",
    "text": "empty(::Type{<:Elements}) -> Tuple{}\n\nGet the empty elements.\n\nThe empty elements is defined to be the empty tuple.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Base.findfirst-Tuple{Hamiltonian.Mathematics.AlgebraOverFields.ID,Hamiltonian.Mathematics.AlgebraOverFields.IdSpace}",
    "page": "Algebra over fields",
    "title": "Base.findfirst",
    "category": "method",
    "text": "findfirst(id::ID,idspace::IdSpace) -> Int\nsearchsortedfirst(idspace::IdSpace,id::ID) -> Int\n\nFind the index of an id in a idspace.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Base.getindex-Tuple{Hamiltonian.Mathematics.AlgebraOverFields.IdSpace,Int64}",
    "page": "Algebra over fields",
    "title": "Base.getindex",
    "category": "method",
    "text": "getindex(idspace::IdSpace,i::Int) -> ID\n\nGet the ith id of a idspace.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Base.getproperty-Tuple{Hamiltonian.Mathematics.AlgebraOverFields.ID,Symbol}",
    "page": "Algebra over fields",
    "title": "Base.getproperty",
    "category": "method",
    "text": "getproperty(cid::ID,name::Symbol)\n\nGet the property of a composite id.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Base.isless-Tuple{Hamiltonian.Mathematics.AlgebraOverFields.ID,Hamiltonian.Mathematics.AlgebraOverFields.ID}",
    "page": "Algebra over fields",
    "title": "Base.isless",
    "category": "method",
    "text": "isless(cid1::ID,cid2::ID) -> Bool\n<(cid1::ID,cid2::ID) -> Bool\n\nCompare two ids and judge whether the first is less than the second.\n\nWe assume that ids with smaller ranks are always less than those with higher ranks. If two ids are of the same rank, the comparison goes just like that between tuples.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Base.propertynames-Union{Tuple{Type{I}}, Tuple{I}, Tuple{Type{I},Bool}} where I<:Hamiltonian.Mathematics.AlgebraOverFields.ID",
    "page": "Algebra over fields",
    "title": "Base.propertynames",
    "category": "method",
    "text": "propertynames(::Type{I},private::Bool=false) where I<:ID -> Tuple\n\nGet the property names of a composite id.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Base.replace-Tuple{Hamiltonian.Mathematics.AlgebraOverFields.Element}",
    "page": "Algebra over fields",
    "title": "Base.replace",
    "category": "method",
    "text": "replace(m::Element;kwargs...) -> typeof(m)\n\nReturn a copy of a concrete Element with some of the field values replaced by the keyword arguments.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Base.show-Tuple{IO,Dict{I,M} where M<:Hamiltonian.Mathematics.AlgebraOverFields.Element where I<:Hamiltonian.Mathematics.AlgebraOverFields.ID}",
    "page": "Algebra over fields",
    "title": "Base.show",
    "category": "method",
    "text": "show(io::IO,ms::Elements)\n\nShow a set of elements.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Base.show-Tuple{IO,Hamiltonian.Mathematics.AlgebraOverFields.ID}",
    "page": "Algebra over fields",
    "title": "Base.show",
    "category": "method",
    "text": "show(io::IO,cid::ID)\n\nShow a composite id.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Base.valtype-Union{Tuple{Type{#s68} where #s68<:(Element{N,V,I} where I<:(ID{#s69} where #s69<:Tuple{Vararg{SimpleID,N}}))}, Tuple{V}, Tuple{N}} where V<:Number where N",
    "page": "Algebra over fields",
    "title": "Base.valtype",
    "category": "method",
    "text": "valtype(::Type{<:Element{N,V}}) where {N,V<:Number}\nvaltype(m::Element)\n\nGet the type of the value of an element.\n\nThe result is also the type of the field over which the algebra is defined.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Base.zero-Tuple{Dict{I,M} where M<:Hamiltonian.Mathematics.AlgebraOverFields.Element where I<:Hamiltonian.Mathematics.AlgebraOverFields.ID}",
    "page": "Algebra over fields",
    "title": "Base.zero",
    "category": "method",
    "text": "zero(ms::Elements) -> typeof(ms)\nzero(::Type{Elements{I,M}}) where {I,M} -> Elements{I,M}\n\nGet a zero set of elements.\n\nA zero set of elements is defined to be the one with no elements.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/AlgebraOverFields.html#Manual-1",
    "page": "Algebra over fields",
    "title": "Manual",
    "category": "section",
    "text": "Modules=[AlgebraOverFields]\nOrder=  [:module,:constant,:type,:macro,:function]"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#",
    "page": "Quantum numbers",
    "title": "Quantum numbers",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Mathematics.QuantumNumbers"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Quantum-numbers-1",
    "page": "Quantum numbers",
    "title": "Quantum numbers",
    "category": "section",
    "text": "Qunatum numbers can be considered as the conserved labels for the bases of a Hilbert space when a quantum system hosts some symmetries. Here we only implement Abelian quantum numbers because non-Abelian ones are far more complicated yet much less used. In practice, quantum numbers can be integers or half integers, therefore, we use real numbers to denote them in this module for simplicity. Independent quantum numbers, such as the particle number and the spin z-component, can coexist at the same time. We use type AbelianNumber to represent the complete set of independent ones for a single basis of a Hilbert space, and type AbelianNumbers to represent the whole quantum numbers for the total bases."
},

{
    "location": "man/Mathematics/QuantumNumbers.html#AbelianNumber-1",
    "page": "Quantum numbers",
    "title": "AbelianNumber",
    "category": "section",
    "text": "The abstract type for the complete set of independent quantum numbers for a single basis.Main features include:function fieldnames: get the names of the quantum numbers\nfunction periods: get the periods of the quantum numbers\narithmetic operations: +, -, *, ^, ⊕, ⊗\nhashable: concrete instances can be used as keys for a dict or a set\niterable: concrete instances are iterable over their values\ncomparable: two concrete instances can be comparedIn particular, AbelianNumber <: HomoNamedVector{Float}, all features supported by HomoNamedVector are also available for HomoNamedVector. See also HomoNamedVector.For convenience, 4 kinds of quantum numbers are predefined in this module, i.e.SQN: for spin z-component reserved systems\nPQN: for particle number reserved systems\nSPQN: for both particle number and spin-z component reserved systems\nZ2QN: for systems with a Z_2 conservation quantum numberUsers who want to define their own Z_N-like quantum numbers must handle the periodicities in the construction function, otherwise, wrong results will be get when arithmetic operations, such as + or -, are involved. It is highly recommended to use the macro @quantumnumber to define your own concrete AbelianNumbers."
},

{
    "location": "man/Mathematics/QuantumNumbers.html#AbelianNumbers-1",
    "page": "Quantum numbers",
    "title": "AbelianNumbers",
    "category": "section",
    "text": "The whole quantum numbers for the total bases.By design, a AbelianNumbers{QN} has one type parameter:QN<:AbelianNumber: the type of the quantum numbers contained in itAnd 3 attributes:form::Char: Its form, whose value must be one of the followings\n\'G\': the general form, which has no restriction for its contents\n\'U\': the unitary form, which requires no duplicates in its contents\n\'C\': the canonical form, which requires both no duplicates and accending-order in its contents\nUsually, G-formed and U-formed AbelianNumberses can be transformed to the corresponding C-formed ones by the sort function.\ncontents::Vector{QN}: The quantum numbers contained in it. To achieve high efficiency, it is required to be an homogenous array of a certain kind of concrete AbelianNumber.\nindptr::Vector{Int}: The indptr of the quantum numbers contained in it, which is similar to the colptr attribute of a CSC sparse matrix and records the compression info of its contents.Main features include:function eltype: get the concrete type of the quantum numbers it contains\nindex access: get the contents directly by the getindex function\narithmetic operations: +, -, *, ^, ⊗, ⊕\niterable: various iteration supports, including functions such as iterate, keys, values and pairsFor a complete summation of its features, please refer to the manual.For convenience, 5 functions are predefined to generate the AbelianNumbers of common physical systems, i.e.SQNS: a signle spin\nPQNS: a single-particle state with at most N identical particles\nSzPQNS: a single-paritcle state with at most one particle whose spin-z component is Sz\nSPQNS: a single site with internal degrees of freedom that can be ascribed to a spin\nZ2QNS: any Z_2 Hilbert space"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Mathematics.QuantumNumbers.qnsbruteforce",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Mathematics.QuantumNumbers.qnsbruteforce",
    "category": "constant",
    "text": "qnsbruteforce\n\nIndicate that decompose uses the brute force method.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Mathematics.QuantumNumbers.qnscompression",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Mathematics.QuantumNumbers.qnscompression",
    "category": "constant",
    "text": "qnscompression\n\nIndicate that findall and permute use the compressed contents.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Mathematics.QuantumNumbers.qnscontents",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Mathematics.QuantumNumbers.qnscontents",
    "category": "constant",
    "text": "qnscontents\n\nIndicate that expand uses the compressed/expanded contents.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Mathematics.QuantumNumbers.qnscounts",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Mathematics.QuantumNumbers.qnscounts",
    "category": "constant",
    "text": "qnscounts\n\nIndicate that methods with AbelianNumbers use the count number of the compressed contents.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Mathematics.QuantumNumbers.qnsexpansion",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Mathematics.QuantumNumbers.qnsexpansion",
    "category": "constant",
    "text": "qnsexpansion\n\nIndicate that findall and permute use the expanded contents.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Mathematics.QuantumNumbers.qnsindices",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Mathematics.QuantumNumbers.qnsindices",
    "category": "constant",
    "text": "qnsindices\n\nIndicate that expand uses the indices of the compressed/expanded contents.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Mathematics.QuantumNumbers.qnsindptr",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Mathematics.QuantumNumbers.qnsindptr",
    "category": "constant",
    "text": "qnsindptr\n\nIndicate that methods with AbelianNumbers use the index pointer of the compressed contents.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Mathematics.QuantumNumbers.qnsmontecarlo",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Mathematics.QuantumNumbers.qnsmontecarlo",
    "category": "constant",
    "text": "qnsmontecarlo\n\nIndicate that decompose uses the Monte Carlo method.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Mathematics.QuantumNumbers.AbelianNumber",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Mathematics.QuantumNumbers.AbelianNumber",
    "category": "type",
    "text": "Abstract type for all concrete quantum numbers for a single basis.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Mathematics.QuantumNumbers.AbelianNumbers",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Mathematics.QuantumNumbers.AbelianNumbers",
    "category": "type",
    "text": "AbelianNumbers(form::Char,contents::Vector{<:AbelianNumber},counts::Vector{Int},::QnsCounts)\nAbelianNumbers(form::Char,contents::Vector{<:AbelianNumber},indptr::Vector{Int},::QnsIndptr)\n\nThe whole quantum numbers of the total bases of a Hilbert space.\n\nThe default constructors construct a AbelianNumbers from a vector of concrete quantum numbers and an vector containing their counts or indptr.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Mathematics.QuantumNumbers.AbelianNumbers",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Mathematics.QuantumNumbers.AbelianNumbers",
    "category": "type",
    "text": "AbelianNumbers(qn::AbelianNumber,count::Int=1)\n\nConstruct a AbelianNumbers with one unique quantum number which occurs count times.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Mathematics.QuantumNumbers.AbelianNumbers-Tuple{OrderedCollections.OrderedDict{#s15,Int64} where #s15<:Hamiltonian.Mathematics.QuantumNumbers.AbelianNumber}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Mathematics.QuantumNumbers.AbelianNumbers",
    "category": "method",
    "text": "AbelianNumbers(od::OrderedDict{<:AbelianNumber,Int})\n\nConstruct a AbelianNumbers from an ordered dict containing concrete quantum numbers and their counts.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Mathematics.QuantumNumbers.AbelianNumbers-Tuple{OrderedCollections.OrderedDict{#s15,UnitRange{Int64}} where #s15<:Hamiltonian.Mathematics.QuantumNumbers.AbelianNumber}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Mathematics.QuantumNumbers.AbelianNumbers",
    "category": "method",
    "text": "AbelianNumbers(od::OrderedDict{<:AbelianNumber,UnitRange{Int}})\n\nConstruct a AbelianNumbers from an ordered dict containing concrete quantum numbers and their slices.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Mathematics.QuantumNumbers.PQN",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Mathematics.QuantumNumbers.PQN",
    "category": "type",
    "text": "PQN(N::Real)\n\nThe concrete AbelianNumber of a quantum system with particle number N conserved.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Mathematics.QuantumNumbers.SPQN",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Mathematics.QuantumNumbers.SPQN",
    "category": "type",
    "text": "SPQN(N::Real,Sz::Real)\n\nThe concrete AbelianNumber of a quantum system with both particle number N and spin z-component Sz conserved.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Mathematics.QuantumNumbers.SQN",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Mathematics.QuantumNumbers.SQN",
    "category": "type",
    "text": "SQN(Sz::Real)\n\nThe concrete AbelianNumber of a quantum system with spin z-component Sz conserved.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Mathematics.QuantumNumbers.Z2QN",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Mathematics.QuantumNumbers.Z2QN",
    "category": "type",
    "text": "Z2QN(N::Real)\n\nThe concrete AbelianNumber of a quantum system with a Z₂-like conserved quantity.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Mathematics.QuantumNumbers.@quantumnumber-Tuple{Any,Any,Any}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Mathematics.QuantumNumbers.@quantumnumber",
    "category": "macro",
    "text": "@quantumnumber typename fieldnames fieldperiods\n\nConstruct a concrete AbelianNumber with the type name being typename, fieldnames specified by fieldnames and periods specified by fieldperiods.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Mathematics.QuantumNumbers.PQNS-Tuple{Real}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Mathematics.QuantumNumbers.PQNS",
    "category": "method",
    "text": "PQNS(N::Real) -> AbelianNumbers{PQN}\n\nConstruct the AbelianNumbers of the Hilbert space of a single-particle state with at most N identical particles.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Mathematics.QuantumNumbers.SPQNS-Tuple{Real}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Mathematics.QuantumNumbers.SPQNS",
    "category": "method",
    "text": "SPQNS(S::Real) -> AbelianNumbers{SPQN}\n\nConstruct the AbelianNumbers of the Hilbert space of a single site with internal degrees of freedom that can be ascribed to a spin S.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Mathematics.QuantumNumbers.SQNS-Tuple{Real}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Mathematics.QuantumNumbers.SQNS",
    "category": "method",
    "text": "SQNS(S::Real) -> AbelianNumbers{SQN}\n\nConstruct the AbelianNumbers of the Hilbert space of a signle spin S.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Mathematics.QuantumNumbers.SzPQNS-Tuple{Real}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Mathematics.QuantumNumbers.SzPQNS",
    "category": "method",
    "text": "SzPQNS(Sz::Real) -> AbelianNumbers{SPQN}\n\nConstruct the AbelianNumbers of the Hilbert space of a single-paritcle state with at most one particle whose spin-z component is Sz.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Mathematics.QuantumNumbers.Z2QNS-Tuple{}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Mathematics.QuantumNumbers.Z2QNS",
    "category": "method",
    "text": "Z2QNS() -> AbelianNumbers{Z2QN}\n\nConstruct the AbelianNumbers of a Z_2 Hilbert space.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Mathematics.QuantumNumbers.decompose-Union{Tuple{QN}, Tuple{N}, Tuple{Tuple{Vararg{AbelianNumbers{QN},N}},QN,Tuple{Vararg{Int64,N}},QnsBruteForce}} where QN<:Hamiltonian.Mathematics.QuantumNumbers.AbelianNumber where N",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Mathematics.QuantumNumbers.decompose",
    "category": "method",
    "text": "decompose(qnses::NTuple{N,AbelianNumbers{QN}},target::QN,signs::NTuple{N,Int},::QnsBruteForce;nmax::Int=20) where {N,QN<:AbelianNumber} -> Vector{NTuple{N,Int}}\ndecompose(qnses::NTuple{N,AbelianNumbers{QN}},target::QN,signs::NTuple{N,Int},::QnsMonteCarlo;nmax::Int=20) where {N,QN<:AbelianNumber} -> Vector{NTuple{N,Int}}\n\nFind a couple of decompositions of target with respect to qnses.\n\nnote: Note\nA tuple of integers (i₁,i₂,...) is called a decomposition of a given target with respect to the given qnses if and only if they satisfy the \"decomposition rule\":sum_textj textsignstextjtimestextqnsestextjtexti_textj==texttargetThis equation is in fact a kind of a set of restricted linear Diophantine equations. Indeed, our quantum numbers are always discrete Abelian ones and all instances of a concrete AbelianNumber forms a module over the ring of integers. Therefore, each quantum number can be represented as a integral multiple of the unit element of the Abelian module, which results in the final reduction of the above equation to a set of linear Diophantine equations. Then finding a decomposition is equivalent to find a solution of the reduced linear Diophantine equations, with the restriction that the quantum numbers constructed from the solution should be in the corresponding qnses. Here we provide two methods to find such decompositions, one is by brute force (qnsbruteforce case), and the other is by Monte Carlo simultatioins (qnsmontecarlo case).\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Mathematics.QuantumNumbers.regularize!-Union{Tuple{QN}, Tuple{Type{QN},AbstractArray{Float64,1}}} where QN<:Hamiltonian.Mathematics.QuantumNumbers.AbelianNumber",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Mathematics.QuantumNumbers.regularize!",
    "category": "method",
    "text": "regularize!(::Type{QN},array::AbstractVector{Float}) where QN<:AbelianNumber -> typeof(array)\nregularize!(::Type{QN},array::AbstractMatrix{Float}) where QN<:AbelianNumber -> typeof(array)\n\nRegularize the elements of an array in place so that it can represent quantum numbers.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Mathematics.QuantumNumbers.regularize-Union{Tuple{QN}, Tuple{Type{QN},Union{AbstractArray{Float64,1}, AbstractArray{Float64,2}}}} where QN<:Hamiltonian.Mathematics.QuantumNumbers.AbelianNumber",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Mathematics.QuantumNumbers.regularize",
    "category": "method",
    "text": "regularize(::Type{QN},array::Union{AbstractVector{Float},AbstractMatrix{Float}}) where {QN<:AbelianNumber} -> typeof(array)\n\nRegularize the elements of an array and return a copy that can represent quantum numbers.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Mathematics.QuantumNumbers.toordereddict-Tuple{Hamiltonian.Mathematics.QuantumNumbers.AbelianNumbers,Hamiltonian.Mathematics.QuantumNumbers.QnsIndptr}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Mathematics.QuantumNumbers.toordereddict",
    "category": "method",
    "text": "toordereddict(qns::AbelianNumbers,::QnsIndptr) -> OrderedDict{qns|>eltype,UnitRange{Int}}\ntoordereddict(qns::AbelianNumbers,::QnsCounts) -> OrderedDict{qns|>eltype,Int}\n\nConvert a AbelianNumbers to an ordered dict.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Mathematics.QuantumNumbers.ukron-Union{Tuple{Vararg{AbelianNumbers{QN},N}}, Tuple{QN}, Tuple{N}} where QN<:Hamiltonian.Mathematics.QuantumNumbers.AbelianNumber where N",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Mathematics.QuantumNumbers.ukron",
    "category": "method",
    "text": "ukron(qnses::Vararg{AbelianNumbers{QN},N};signs::NTuple{N,Int}=ntuple(i->1,N)) where {N,QN<:AbelianNumber} -> AbelianNumbers{QN},Dict{QN,Dict{NTuple{N,QN},UnitRange{Int}}}\n\nUnitary Kronecker product of several AbelianNumberses. The product result as well as the records of the product will be returned.\n\nnote: Note\nAll input AbelianNumbers must be \'U\' formed or \'C\' formed.\nSince duplicate quantum number are not allowed in \'U\' formed and \'C\' formed AbelianNumberses, in general, there exists a merge process of duplicate quantum numbers in the product result. Therefore, records are needed to keep track of this process, which will be returned along with the product result. The records are stored in a Dict{QN,Dict{NTuple{N,QN},UnitRange{Int}}} typed dict, in which, for each unduplicate quantum number qn in the product result, there exist a record Dict((qn₁,qn₂,...)=>start:stop,...) telling what quantum numbers (qn₁,qn₂,...) a mereged duplicate qn comes from and what slice start:stop this merged duplicate corresponds.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Prerequisites.Interfaces.:⊕-Tuple{Vararg{Hamiltonian.Mathematics.QuantumNumbers.AbelianNumber,N} where N}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Prerequisites.Interfaces.:⊕",
    "category": "method",
    "text": "⊕(qns::AbelianNumber...) -> AbelianNumbers{qns|>eltype}\n⊕(qnses::AbelianNumbers...) -> qnses|>eltype\n\nGet the direct sum of some AbelianNumbers or AbelianNumberses.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Prerequisites.Interfaces.:⊗-Tuple{Vararg{Hamiltonian.Mathematics.QuantumNumbers.AbelianNumber,N} where N}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Prerequisites.Interfaces.:⊗",
    "category": "method",
    "text": "⊗(qns::AbelianNumber...) -> eltype(qns)\n⊗(qnses::AbelianNumbers...) -> eltype(qnses)\n\nGet the direct product of some AbelianNumbers or AbelianNumberses.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Prerequisites.Interfaces.dimension-Tuple{Hamiltonian.Mathematics.QuantumNumbers.AbelianNumbers}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Prerequisites.Interfaces.dimension",
    "category": "method",
    "text": "dimension(qns::AbelianNumbers) -> Int\n\nThe dimension of the Hilbert space a AbelianNumbers represents.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Prerequisites.Interfaces.dimension-Tuple{Type{#s68} where #s68<:Hamiltonian.Mathematics.QuantumNumbers.AbelianNumber}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Prerequisites.Interfaces.dimension",
    "category": "method",
    "text": "dimension(::Type{<:AbelianNumber}) -> Int\ndimension(::AbelianNumber) -> Int\n\nThe dimension of the Hilbert space a AbelianNumber represents. Apparently, this is always 1.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Prerequisites.Interfaces.expand-Tuple{Hamiltonian.Mathematics.QuantumNumbers.AbelianNumbers,Hamiltonian.Mathematics.QuantumNumbers.QnsContents}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Prerequisites.Interfaces.expand",
    "category": "method",
    "text": "expand(qns::AbelianNumbers,::QnsContents) -> Vector{qns|>eltype}\nexpand(qns::AbelianNumbers,::QnsIndices) -> Vector{Int}\n\nExpand the contents (qnscontents case) or indices (qnsindices case) of a AbelianNumbers to the uncompressed form.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Hamiltonian.Prerequisites.Interfaces.permute-Tuple{Hamiltonian.Mathematics.QuantumNumbers.AbelianNumbers,Array{Int64,1},Hamiltonian.Mathematics.QuantumNumbers.QnsCompression}",
    "page": "Quantum numbers",
    "title": "Hamiltonian.Prerequisites.Interfaces.permute",
    "category": "method",
    "text": "permute(qns::AbelianNumbers,permutation::Vector{Int},::QnsCompression) -> AbelianNumbers\npermute(qns::AbelianNumbers,permutation::Vector{Int},::QnsExpansion) -> AbelianNumbers\n\nReorder the quantum numbers contained in a AbelianNumbers with a permutation and return the new one.\n\nFor qnscompression case, the permutation is for the compressed contents of the original AbelianNumbers while for qnsexpansion case, the permutation is for the expanded contents of the original AbelianNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Base.:*-Tuple{Hamiltonian.Mathematics.QuantumNumbers.AbelianNumber,Integer}",
    "page": "Quantum numbers",
    "title": "Base.:*",
    "category": "method",
    "text": "*(qn::AbelianNumber,factor::Integer) -> typeof(qn)\n*(factor::Integer,qn::AbelianNumber) -> typeof(qn)\n*(qns::AbelianNumbers,factor::Integer) -> AbelianNumbers\n*(factor::Integer,qns::AbelianNumbers) -> AbelianNumbers\n\nOverloaded * operator for the multiplication between an integer and a AbelianNumber or a AbelianNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Base.:+-Tuple{Hamiltonian.Mathematics.QuantumNumbers.AbelianNumber}",
    "page": "Quantum numbers",
    "title": "Base.:+",
    "category": "method",
    "text": "+(qn::AbelianNumber) -> typeof(qn)\n+(qn::QN,qns::QN...) where QN<:AbelianNumber -> QN\n+(qns::AbelianNumbers) -> AbelianNumbers\n+(qn::QN,qns::AbelianNumbers{QN}) where QN<:AbelianNumber -> AbelianNumbers{QN}\n+(qns::AbelianNumbers{QN},qn::QN) where QN<:AbelianNumber -> AbelianNumbers{QN}\n\nOverloaded + operator for AbelianNumber and AbelianNumbers.\n\nnote: Note\nThe addition between a AbelianNumbers and an AbelianNumber is just a global shift of the contents of the AbelianNumbers by the AbelianNumber, therefore, the result is a AbelianNumbers.\n+ cannot be used between two AbelianNumbers because the result is ambiguous. Instead, use ⊕ for direct sum and ⊗ for direct product.\nTo ensure type stability, two AbelianNumber can be added together if and only if they are of the same type.\nSimilarly, a AbelianNumber and a AbelianNumbers can be added together if and only if the former\'s type is the same with the latter\'s eltype.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Base.:--Tuple{Hamiltonian.Mathematics.QuantumNumbers.AbelianNumber}",
    "page": "Quantum numbers",
    "title": "Base.:-",
    "category": "method",
    "text": "-(qn::AbelianNumber) -> typeof(qn)\n-(qn1::QN,qn2::QN) where QN<:AbelianNumber -> QN\n-(qns::AbelianNumbers) -> AbelianNumbers\n-(qn::QN,qns::AbelianNumbers{QN}) where QN<:AbelianNumber -> AbelianNumbers{QN}\n-(qns::AbelianNumbers{QN},qn::QN) where QN<:AbelianNumber -> AbelianNumbers{QN}\n\nOverloaded - operator for AbelianNumber and AbelianNumbers.\n\nnote: Note\nThe subtraction between a AbelianNumbers and a AbelianNumber is just a global shift of the contents of the AbelianNumbers by the AbelianNumber, therefore, the result is a AbelianNumbers.\n- cannot be used between two AbelianNumbers because the result is ambiguous. Instead, use ⊕ with signs for direct sum and ⊗ with signs for direct product.\nTo ensure type stability, a AbelianNumber can be subtracted by another AbelianNumber if and only if they are of the same type.\nSimilarly, a AbelianNumber can be subtracted by a AbelianNumbers or vice versa if and only if the former\'s type is the same with the latter\'s eltype.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Base.:==-Tuple{Hamiltonian.Mathematics.QuantumNumbers.AbelianNumbers,Hamiltonian.Mathematics.QuantumNumbers.AbelianNumbers}",
    "page": "Quantum numbers",
    "title": "Base.:==",
    "category": "method",
    "text": "==(qns1::AbelianNumbers,qns2::AbelianNumbers) -> Bool\nisequal(qns1::AbelianNumbers,qns2::AbelianNumbers) -> Bool\n\nOverloaded equivalent operator. Two AbelianNumberses are equal to each other if and only if both their contentses and indptrs are elementwise equal to each other.\n\nnote: Note\nIt is not necessary for two AbelianNumberses to have the same eltype nor the same form to be equal to each other.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Base.:^-Tuple{Hamiltonian.Mathematics.QuantumNumbers.AbelianNumber,Integer}",
    "page": "Quantum numbers",
    "title": "Base.:^",
    "category": "method",
    "text": "^(qn::AbelianNumber,factor::Integer) -> typeof(qn)\n^(qns::AbelianNumbers,factor::Integer) -> AbelianNumbers\n\nOverloaded ^ operator for AbelianNumber and AbelianNumbers. This operation translates into the direct product of factor copies of qn or qns.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Base.eltype-Union{Tuple{Type{#s68} where #s68<:AbelianNumbers{QN}}, Tuple{QN}} where QN",
    "page": "Quantum numbers",
    "title": "Base.eltype",
    "category": "method",
    "text": "eltype(::Type{<:AbelianNumbers{QN}}) where QN\neltype(qns::AbelianNumbers)\n\nGet the type of the concrete AbelianNumber contained in a AbelianNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Base.filter-Union{Tuple{QN}, Tuple{QN,AbelianNumbers{QN}}} where QN<:Hamiltonian.Mathematics.QuantumNumbers.AbelianNumber",
    "page": "Quantum numbers",
    "title": "Base.filter",
    "category": "method",
    "text": "filter(target::QN,qns::AbelianNumbers{QN}) where QN<:AbelianNumber -> AbelianNumbers{QN}\nfilter(targets::NTuple{N,QN},qns::AbelianNumbers{QN}) where {N,QN<:AbelianNumber} -> AbelianNumbers{QN}\n\nFind a subset of a AbelianNumbers by picking out the quantum numbers in targets.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Base.findall-Union{Tuple{QN}, Tuple{QN,AbelianNumbers{QN},QnsCompression}} where QN<:Hamiltonian.Mathematics.QuantumNumbers.AbelianNumber",
    "page": "Quantum numbers",
    "title": "Base.findall",
    "category": "method",
    "text": "findall(target::QN,qns::AbelianNumbers{QN},::QnsCompression) where QN<:AbelianNumber -> Vector{Int}\nfindall(target::QN,qns::AbelianNumbers{QN},::QnsExpansion) where QN<:AbelianNumber -> Vector{Int}\nfindall(targets::NTuple{N,QN},qns::AbelianNumbers{QN},::QnsCompression) where {N,QN<:AbelianNumber} -> Vector{Int}\nfindall(targets::NTuple{N,QN},qns::AbelianNumbers{QN},::QnsExpansion) where {N,QN<:AbelianNumber} -> Vector{Int}\n\nFind all the indices of the target quantum numbers in the contents (qnscompression case) or the expansion (qnsexpansion case) of a AbelianNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Base.getindex-Tuple{Hamiltonian.Mathematics.QuantumNumbers.AbelianNumbers,Int64}",
    "page": "Quantum numbers",
    "title": "Base.getindex",
    "category": "method",
    "text": "getindex(qns::AbelianNumbers,index::Int) -> eltype(qns)\ngetindex(qns::AbelianNumbers,slice::UnitRange{Int}) -> AbelianNumbers\ngetindex(qns::AbelianNumbers,indices::Vector{Int}) -> AbelianNumbers\n\nOverloaded [] operator.\n\nnote: Note\nFor a AbelianNumbers, all these getindex functions act on its contents, i.e. its compressed data, but not on its expansion, i.e. the uncompressed data. This definition is consistent with the length function.\nWhen the index is an integer, the result is a AbelianNumber, while when the index is a unit range or a vector of intgers, the result is a AbelianNumbers. The logic is quite reasonable because such behaviors are much alike to those of a vector container.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Base.iterate",
    "page": "Quantum numbers",
    "title": "Base.iterate",
    "category": "function",
    "text": "iterate(qns::AbelianNumbers,state::Int=1)\niterate(rv::Iterators.Reverse{<:AbelianNumbers},state::Int=length(rv.itr,false))\n\nIterate or reversely iterate over the concrete AbelianNumbers contained in a AbelianNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Base.keys-Tuple{Hamiltonian.Mathematics.QuantumNumbers.AbelianNumbers}",
    "page": "Quantum numbers",
    "title": "Base.keys",
    "category": "method",
    "text": "keys(qns::AbelianNumbers) -> Vector{qns|>eltype}\n\nIterate over the concrete AbelianNumbers contained in a AbelianNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Base.kron-Union{Tuple{QN}, Tuple{Type{QN},AbelianNumber,AbelianNumber}} where QN<:Hamiltonian.Mathematics.QuantumNumbers.AbelianNumber",
    "page": "Quantum numbers",
    "title": "Base.kron",
    "category": "method",
    "text": "kron(::Type{QN},qn1::AbelianNumber,qn2::AbelianNumber) where QN<:AbelianNumber -> QN\nkron(qns::Vararg{<:AbelianNumber,N};signs::NTuple{N,Int}=ntuple(i->1,N)) where N -> eltype(qns)\nkron(qnses::Vararg{AbelianNumbers{QN},N};signs::NTuple{N,Int}=ntuple(i->1,N)) where {N,QN<:AbelianNumber} -> AbelianNumbers{QN}\n\nGet the direct product of some AbelianNumbers or AbelianNumberses.\n\nnote: Note\nPhysically, the direct product of a couple of AbelianNumbers or AbelianNumberses are defined by the direct product of the bases of the Hilbert spaces they represent. Therefore, AbelianNumbers with differenct types or AbelianNumberses with differenct eltypes are allowed to be direct producted in principle. However, for simplicity, we only implement a method which handle the situation of two AbelianNumbers with differenct types. The type of the result should be provided as the first parameter. Note that in this situation, the fieldnames and periods of the result type must be exactly equal to the flattened fieldnames and periods of the two input AbelianNumbers, which means, even the order of the input AbelianNumbers matters.\nApparently, the dimension of the result equals the product of those of the inputs. Therefore, the direct product of AbelianNumbers is also a AbelianNumber since its dimension is still one.\nFor other situations except the one mentioned in Note.1, the input AbelianNumbers or AbelianNumberses must be homogenous. Meanwhile, signs can also be provided for these situations. Note that each quantum number in the contents of the result is obtained by a summation of the corresponding quanum numbers out of the inputs with the correct signs. This is a direct observation of the Abelian nature of our quantum numbers.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Base.length-Tuple{Hamiltonian.Mathematics.QuantumNumbers.AbelianNumbers}",
    "page": "Quantum numbers",
    "title": "Base.length",
    "category": "method",
    "text": "length(qns::AbelianNumbers) -> Int\n\nGet the number of unduplicate qunatum numbers in the AbelianNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Base.pairs-Tuple{Hamiltonian.Mathematics.QuantumNumbers.AbelianNumbers,Union{QnsCounts, QnsIndptr}}",
    "page": "Quantum numbers",
    "title": "Base.pairs",
    "category": "method",
    "text": "pairs(qns::AbelianNumbers,choice::Union{QnsIndptr,QnsCounts})\n\nIterate over the AbelianNumber=>slice or AbelianNumber=>count pairs.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Base.show-Tuple{IO,Hamiltonian.Mathematics.QuantumNumbers.AbelianNumbers}",
    "page": "Quantum numbers",
    "title": "Base.show",
    "category": "method",
    "text": "show(io::IO,qns::AbelianNumbers)\n\nShow a AbelianNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Base.sort-Tuple{Hamiltonian.Mathematics.QuantumNumbers.AbelianNumbers}",
    "page": "Quantum numbers",
    "title": "Base.sort",
    "category": "method",
    "text": "sort(qns::AbelianNumbers) -> AbelianNumbers,Vector{Int}\n\nSort the quantum numbers of a AbelianNumber, return the sorted AbelianNumber and the permutation array that sorts the expansion of the original AbelianNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Base.string-Tuple{Hamiltonian.Mathematics.QuantumNumbers.AbelianNumbers}",
    "page": "Quantum numbers",
    "title": "Base.string",
    "category": "method",
    "text": "string(qns::AbelianNumbers) -> String\n\nConvert a AbelianNumbers to string.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Base.union-Union{Tuple{Vararg{AbelianNumber,N}}, Tuple{N}} where N",
    "page": "Quantum numbers",
    "title": "Base.union",
    "category": "method",
    "text": "union(qns::Vararg{<:AbelianNumber,N};signs::NTuple{N,Int}=ntuple(i->1,N)) where N -> AbelianNumbers\nunion(qnses::Vararg{AbelianNumbers{QN},N};signs::NTuple{N,Int}=ntuple(i->1,N)) where {N,QN<:AbelianNumber} -> AbelianNumbers{QN}\n\nGet the direct sum of some AbelianNumbers or AbelianNumberses.\n\nnote: Note\nPhysically, the direct sum of a couple of AbelianNumbers or AbelianNumberses is defined by the direct sum of the bases of the Hilbert spaces they represent. Therefore, the input AbelianNumbers or AbelianNumberses must be homogenous. Inhomogenous \'AbelianNumber\'s must be direct producted first to ensure homogenity before the direct sum.\nApparently, the dimension of the result equals the summation of those of the inputs, which means, even for AbelianNumbers, the result will be naturally a AbelianNumbers because the dimension of the result is larger than 1.\nSigns of AbelianNumbers or AbelianNumberses can be provided when getting their direct sums.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#Base.values-Tuple{Hamiltonian.Mathematics.QuantumNumbers.AbelianNumbers,Hamiltonian.Mathematics.QuantumNumbers.QnsIndptr}",
    "page": "Quantum numbers",
    "title": "Base.values",
    "category": "method",
    "text": "values(qns::AbelianNumbers,::QnsIndptr)\nvalues(qns::AbelianNumbers,::QnsCounts)\n\nIterate over the slices/counts of the AbelianNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Mathematics/QuantumNumbers.html#qnmanual-1",
    "page": "Quantum numbers",
    "title": "Manual",
    "category": "section",
    "text": "Modules=[QuantumNumbers]\nOrder=  [:module,:constant,:type,:macro,:function]"
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
    "text": "Essentials of the Hamiltonian package, which defines all the imported constants, types and functions when using import Hamiltonian or using Hamiltonian. Note that this submodule depends on the Prerequisites and Mathematics submodules although the variables in them are not exported to the scope of Hamiltonian by default.Pages=  [\n        \"Spatials.md\",\n        \"DegreesOfFreedom.md\",\n        \"Terms.md\",\n        \"FockPackage.md\",\n        \"SpinPackage.md\",\n        ]\nDepth=2"
},

{
    "location": "man/Essentials/Spatials.html#",
    "page": "Spatials",
    "title": "Spatials",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Essentials.Spatialspush!(LOAD_PATH,\"../../../../src/\")\nusing Hamiltonian"
},

{
    "location": "man/Essentials/Spatials.html#Spatials-1",
    "page": "Spatials",
    "title": "Spatials",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/Spatials.html#Utilities-1",
    "page": "Spatials",
    "title": "Utilities",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/Spatials.html#AbstractBond-1",
    "page": "Spatials",
    "title": "AbstractBond",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/Spatials.html#Point-1",
    "page": "Spatials",
    "title": "Point",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/Spatials.html#Bond-1",
    "page": "Spatials",
    "title": "Bond",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/Spatials.html#AbstractLattice-1",
    "page": "Spatials",
    "title": "AbstractLattice",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/Spatials.html#Lattice-1",
    "page": "Spatials",
    "title": "Lattice",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/Spatials.html#SuperLattice-1",
    "page": "Spatials",
    "title": "SuperLattice",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/Spatials.html#Cylinder-1",
    "page": "Spatials",
    "title": "Cylinder",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.acrossbonds",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.acrossbonds",
    "category": "constant",
    "text": "acrossbonds\n\nIndicate that bonds across the unitcell are inquired, which are in fact those across the periodic boundaries.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.insidebonds",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.insidebonds",
    "category": "constant",
    "text": "insidebonds\n\nIndicate that bonds inside the unitcell are inquired, which do not contain those across the periodic boundaries.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.interbonds",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.interbonds",
    "category": "constant",
    "text": "interbonds\n\nIndicate that bonds inter the sublattices are inquired.\n\nnote: Note\nThese bonds do not contain those accorss the periodic boundaries.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.intrabonds",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.intrabonds",
    "category": "constant",
    "text": "intrabonds\n\nIndicate that bonds intra the sublattices are inquired.\n\nnote: Note\nThese bonds do not contain those accorss the periodic boundaries.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.zerothbonds",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.zerothbonds",
    "category": "constant",
    "text": "zerothbonds\n\nIndicate that zeroth bonds, i.e. the points are inquired.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.AbstractBond",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.AbstractBond",
    "category": "type",
    "text": "AbstractBond{R,P<:PID,N}\n\nAbstract bond.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.AbstractLattice",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.AbstractLattice",
    "category": "type",
    "text": "AbstractLattice{P<:PID,N}\n\nAbstract type for all lattices.\n\nIt should have the following attributes\n\nname::String: the name of the lattice\npids::Vector{P}: the pids of the lattice\nrcoords::Matrix{Float}: the rcoords of the lattice\nicoords::Matrix{Float}: the icoords of the lattice\nvectors::Vector{SVector{N,Float}}: the translation vectors of the lattice\nreciprocals::Vector{SVector{N,Float}}: the reciprocals of the lattice\nneighbors::Dict{Int,Float}: the order-distance map of the nearest neighbors of the lattice\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.AbstractLatticeIndex",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.AbstractLatticeIndex",
    "category": "type",
    "text": "AbstractLatticeIndex{I<:Union{<:PID,Int}}\n\nAbstract index type for a lattice.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.Bond",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.Bond",
    "category": "type",
    "text": "Bond(neighbor::Int,spoint::Point,epoint::Point)\n\nA bond in a lattice.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.Cylinder",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.Cylinder",
    "category": "type",
    "text": "Cylinder{P}(    name::String,\n                block::AbstractMatrix{<:Real},\n                translation::AbstractVector{<:Real};\n                vector::Union{AbstractVector{<:Real},Nothing}=nothing,\n                neighbors::Union{Dict{Int,<:Real},Int}=1,\n                ) where P<:PID\n\nCylinder of 1d and quasi 2d lattices.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.Cylinder-Tuple",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.Cylinder",
    "category": "method",
    "text": "(cylinder::Cylinder)(scopes::Any...;coordination::Int=8) -> Lattice\n\nConstruct a lattice from a cylinder with the assigned scopes.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.ICoordIndex",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.ICoordIndex",
    "category": "type",
    "text": "ICoordIndex(index::Union{<:PID,Int})\n\nIndex for getting an icoord of a lattice.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.Lattice",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.Lattice",
    "category": "type",
    "text": "Lattice(    name::String,\n            pids::Vector{<:PID},\n            rcoords::AbstractMatrix{<:Real};\n            icoords::AbstractMatrix{<:Real}=SMatrix{0,0,Float}(),\n            vectors::AbstractVector{<:AbstractVector{<:Real}}=SVector{0,SVector{size(rcoords,1),Float}}(),\n            neighbors::Union{Dict{Int,<:Real},Int}=1,\n            coordination::Int=8\n        )\nLattice(    name::String,\n            points::AbstractVector{<:Point};\n            vectors::AbstractVector{<:AbstractVector{<:Real}}=SVector{0,SVector{points|>eltype|>dimension,Float}}(),\n            neighbors::Union{Dict{Int,<:Real},Int}=1,\n            coordination::Int=8\n            )\nLattice(    name::String,\n            sublattices::AbstractVector{<:AbstractLattice};\n            vectors::AbstractVector{<:AbstractVector{<:Real}}=SVector{0,SVector{sublattices|>eltype|>dimension,Float}}(),\n            neighbors::Union{Dict{Int,<:Real},Int}=1,\n            coordination::Int=8\n            )\n\nSimplest lattice.\n\nA simplest lattice can be construted from its contents, i.e. pids, rcoords and icoords, or from a couple of points, or from a couple of sublattices.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.PID",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.PID",
    "category": "type",
    "text": "PID(scope,site::Int)\nPID(;scope=\"tz\",site::Int=1)\n\nThe id of a point.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.Point",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.Point",
    "category": "type",
    "text": "Point(pid::PID,rcoord::SVector{N,<:Real},icoord::SVector{N,<:Real}) where N\nPoint(pid::PID,rcoord::NTuple{N,<:Real},icoord::NTuple{N,<:Real}=ntuple(i->0.0,N)) where N\nPoint(pid::PID,rcoord::AbstractVector{<:Real},icoord::AbstractVector{<:Real}=zero(SVector{length(rcoord),Float}))\n\nLabeled point.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.PointIndex",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.PointIndex",
    "category": "type",
    "text": "PointIndex(index::Union{<:PID,Int})\n\nIndex for getting a point of a lattice.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.RCoordIndex",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.RCoordIndex",
    "category": "type",
    "text": "RCoordIndex(index::Union{<:PID,Int})\n\nIndex for getting a rcoord of a lattice.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.SuperLattice",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.SuperLattice",
    "category": "type",
    "text": "SuperLattice(   name::String,\n                sublattices::AbstractVector{<:AbstractLattice};\n                vectors::AbstractVector{<:AbstractVector{<:Real}}=SVector{0,SVector{sublattices|>eltype|>dimension,Float}}(),\n                neighbors::Dict{Int,<:Real}=Dict{Int,Float}()\n                )\n\nSuperLattice that is composed of serveral sublattices.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.azimuth-Tuple{AbstractArray{#s221,1} where #s221<:Real}",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.azimuth",
    "category": "method",
    "text": "azimuth(v::AbstractVector{<:Real}) -> Float\n\nGet the azimuth angle in radians of a vector.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.azimuthd-Tuple{AbstractArray{#s221,1} where #s221<:Real}",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.azimuthd",
    "category": "method",
    "text": "azimuthd(v::AbstractVector{<:Real}) -> Float\n\nGet the azimuth angle in degrees of a vector.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.bonds-Tuple{AbstractLattice}",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.bonds",
    "category": "method",
    "text": "bonds(lattice::AbstractLattice) -> Vector{AbstractBond}\nbonds(lattice::AbstractLattice,::ZerothBonds) -> Vector{Point}\nbonds(lattice::AbstractLattice,::InsideBonds) -> Vector{Bond}\nbonds(lattice::AbstractLattice,::AcrossBonds) -> Vector{Bond}\n\nGet the bonds of a lattice.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.bonds-Tuple{SuperLattice,Hamiltonian.Essentials.Spatials.IntraBonds}",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.bonds",
    "category": "method",
    "text": "bonds(lattice::SuperLattice,::IntraBonds) -> Vector{Bond}\nbonds(lattice::SuperLattice,::InterBonds) -> Vector{Bond}\n\nGet the bonds of a superlattice.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.distance-Tuple{AbstractArray{#s210,1} where #s210<:Real,AbstractArray{#s209,1} where #s209<:Real}",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.distance",
    "category": "method",
    "text": "distance(p1::AbstractVector{<:Real},p2::AbstractVector{<:Real}) -> Float\n\nGet the distance between two points.\n\nnote: Note\nCompared to norm(p1-p2), this function avoids the memory allocation for p1-p2, thus is more efficient.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.icoord-Tuple{Bond}",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.icoord",
    "category": "method",
    "text": "icoord(bond::Bond) -> SVector\n\nGet the icoord of the bond.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.interlinks-Tuple{AbstractArray{#s231,2} where #s231<:Real,AbstractArray{#s230,2} where #s230<:Real,Dict{Int64,Float64}}",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.interlinks",
    "category": "method",
    "text": "interlinks(cluster1::AbstractMatrix{<:Real},cluster2::AbstractMatrix{<:Real},neighbors::Dict{Int,Float}) -> Vector{Tuple{Int,Int,Int,SVector{size(cluster1,1),Float}}}\n\nUse kdtree to get the intercluster nearest neighbors.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.intralinks-Union{Tuple{N}, Tuple{AbstractArray{#s226,2} where #s226<:Real,AbstractArray{#s225,1} where #s225<:(AbstractArray{#s224,1} where #s224<:Real),Dict{Int64,Float64}}, Tuple{AbstractArray{#s223,2} where #s223<:Real,AbstractArray{#s222,1} where #s222<:(AbstractArray{#s221,1} where #s221<:Real),Dict{Int64,Float64},Tuple{Vararg{Int64,N}}}} where N",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.intralinks",
    "category": "method",
    "text": "intralinks( cluster::AbstractMatrix{<:Real},\n            vectors::AbstractVector{<:AbstractVector{<:Real}},\n            neighbors::Dict{Int,Float},\n            maxtranslations::NTuple{N,Int}=ntuple(i->length(neighbors),length(vectors))\n            ) where N -> Vector{Tuple{Int,Int,Int,SVector{size(cluster,1),Float}}}\n\nUse kdtree to get the intracluster nearest neighbors.\n\nAs is similar to minimumlengths, when vectors is nonempty, the cluster assumes periodic boundaries. neighbors provides the map between the bond length and the order of nearest neighbors. Note only those with the lengths present in neighbors will be included in the result. maxtranslations determines the maximum number of translations along those directions specified by vectors when the tiled supercluster is construted (See minimumlengths for the explanation of the method for periodic lattices).\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.isintracell-Tuple{Bond}",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.isintracell",
    "category": "method",
    "text": "isintracell(bond::Bond) -> Bool\n\nJudge whether a bond is intra the unit cell of a lattice.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.isintratriangle-Tuple{AbstractArray{#s15,1} where #s15<:Real,AbstractArray{#s14,1} where #s14<:Real,AbstractArray{#s13,1} where #s13<:Real,AbstractArray{#s12,1} where #s12<:Real}",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.isintratriangle",
    "category": "method",
    "text": "isintratriangle(p::AbstractVector{<:Real},\n                p1::AbstractVector{<:Real},\n                p2::AbstractVector{<:Real},\n                p3::AbstractVector{<:Real};\n                vertexes::NTuple{3,Bool}=(true,true,true),\n                edges::NTuple{3,Bool}=(true,true,true),\n                atol::Real=atol,\n                rtol::Real=rtol\n                ) -> Bool\n\nJudge whether a point belongs to the interior of a triangle whose vertexes are p1, \'p2\' and p3 with the give tolerance. vertexes and edges define whether the interior should contain the vertexes or edges, respectively.\n\nnote: Note\nThe vertexes are in the order (p1,p2,p3) and the edges are in the order (p1p2,p2p3,p3p1).\nThe edges do not contain the vertexes.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.isonline-Tuple{AbstractArray{#s15,1} where #s15<:Real,AbstractArray{#s14,1} where #s14<:Real,AbstractArray{#s13,1} where #s13<:Real}",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.isonline",
    "category": "method",
    "text": "isonline(p::AbstractVector{<:Real},p1::AbstractVector{<:Real},p2::AbstractVector{<:Real};ends::Tuple{Bool,Bool}=(true,true),atol::Real=atol,rtol::Real=rtol) -> Bool\n\nJudge whether a point is on a line segment whose end points are p1 and p2 with the given tolerance. ends defines whether the line segment should contain its ends.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.isparallel-Tuple{AbstractArray{#s68,1} where #s68<:Real,AbstractArray{#s67,1} where #s67<:Real}",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.isparallel",
    "category": "method",
    "text": "isparallel(v1::AbstractVector{<:Real},v2::AbstractVector{<:Real};atol::Real=atol,rtol::Real=rtol) -> Int\n\nJudge whether two vectors are parallel to each other with the given tolerance, 0 for not parallel, 1 for parallel and -1 for antiparallel.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.issubordinate-Tuple{AbstractArray{#s225,1} where #s225<:Real,AbstractArray{#s224,1} where #s224<:(AbstractArray{#s223,1} where #s223<:Real)}",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.issubordinate",
    "category": "method",
    "text": "issubordinate(rcoord::AbstractVector{<:Real},vectors::AbstractVector{<:AbstractVector{<:Real}};atol::Real=atol,rtol::Real=rtol) -> Bool\n\nJudge whether a coordinate belongs to a lattice defined by vectors with the given tolerance.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.minimumlengths",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.minimumlengths",
    "category": "function",
    "text": "minimumlengths(cluster::AbstractMatrix{<:Real},vectors::AbstractVector{<:AbstractVector{<:Real}},nneighbor::Int=1;coordination::Int=8) -> Vector{Float}\n\nUse kdtree to search the lowest several minimum bond lengths within a lattice translated by a cluster.\n\nWhen the translation vectors are not empty, the lattice will be considered periodic in the corresponding directions. Otherwise the lattice will be open in all directions. To search for the bonds accorss the periodic boundaries, the cluster will be pretranslated to become a supercluster, which has open boundaries but is large enough to contain all the nearest neighbors within the required order. The coordination parameter sets the average number of each order of nearest neighbors. If it is to small, larger bond lengths may not be searched, and the result will contain Inf. This is a sign that you may need a larger coordination. Another situation that Inf appears in the result occurs when the minimum lengths are searched in open lattices. Indeed, the cluster may be too small so that the required order just goes beyond it. In this case the warning message can be safely ignored.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.neighbor-Tuple{Bond}",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.neighbor",
    "category": "method",
    "text": "neighbor(bond::Bond) -> Int\n\nGet the neighbor of a bond.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.neighbor-Tuple{Point}",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.neighbor",
    "category": "method",
    "text": "neighbor(::Point) -> 0\nneighbor(::Type{<:Point}) -> 0\n\nGet the neighbor of a point.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.nneighbor-Tuple{AbstractLattice}",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.nneighbor",
    "category": "method",
    "text": "nneighbor(lattice::AbstractLattice) -> Int\n\nGet the highest order of nearest neighbors.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.pidtype-Tuple{AbstractBond}",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.pidtype",
    "category": "method",
    "text": "pidtype(bond::AbstractBond)\npidtype(::Type{<:AbstractBond{R,P}}) where {R,P<:PID}\n\nGet the pid type of a concrete bond.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.polar-Tuple{AbstractArray{#s221,1} where #s221<:Real}",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.polar",
    "category": "method",
    "text": "polar(v::AbstractVector{<:Real}) -> Float\n\nGet the polar angle in radians of a vector.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.polard-Tuple{AbstractArray{#s221,1} where #s221<:Real}",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.polard",
    "category": "method",
    "text": "polard(v::AbstractVector{<:Real}) -> Float\n\nGet the polar angle in degrees of a vector.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.rcoord-Tuple{Bond}",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.rcoord",
    "category": "method",
    "text": "rcoord(bond::Bond) -> SVector\n\nGet the rcoord of the bond.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.reciprocals-Tuple{AbstractArray{#s235,1} where #s235<:(AbstractArray{#s234,1} where #s234<:Real)}",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.reciprocals",
    "category": "method",
    "text": "reciprocals(vectors::AbstractVector{AbstractVector{<:Real}}) -> Vector{Vector{Float}}\n\nGet the reciprocals dual to the input vectors.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.rotate-Tuple{AbstractArray{#s224,2} where #s224<:Real,Real}",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.rotate",
    "category": "method",
    "text": "rotate(cluster::AbstractMatrix{<:Real},angle::Real;axis::Tuple{Union{AbstractVector{<:Real},Nothing},Tuple{<:Real,<:Real}}=(nothing,(0,0))) -> Matrix{Float}\n\nGet a rotated cluster of the original one by a certain angle around an axis.\n\nThe axis is determined by a point it gets through (nothing can be used to denote the origin), and its polar as well as azimuth angles in radians. The default axis is the z axis.\n\nnote: Note\nThe result is given by the Rodrigues\' rotation formula.\nOnly 2 and 3 dimensional vectors can be rotated.\nWhen the input vectors are 2 dimensional, both the polar and azimuth of the axis must be 0.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.tile-Union{Tuple{M}, Tuple{N}, Tuple{AbstractArray{#s226,2} where #s226<:Real,AbstractArray{#s225,1} where #s225<:(AbstractArray{#s224,1} where #s224<:Real)}, Tuple{AbstractArray{#s223,2} where #s223<:Real,AbstractArray{#s222,1} where #s222<:(AbstractArray{#s221,1} where #s221<:Real),Tuple{Vararg{Tuple{Vararg{#s217,N}} where #s217<:Real,M}}}} where M where N",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.tile",
    "category": "method",
    "text": "tile(cluster::AbstractMatrix{<:Real},vectors::AbstractVector{<:AbstractVector{<:Real}},translations::NTuple{M,NTuple{N,<:Real}}=()) where {N,M} -> Matrix{Float}\n\nTile a supercluster by translations of the input cluster.\n\nBasically, the final supercluster is composed of several parts, each of which is a translation of the original cluster, with the translation vectors specified by vectors and each set of the translation indices contained in translations. When translation vectors are empty, a copy of the original cluster will be returned.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.translate-Tuple{AbstractArray{#s235,2} where #s235<:Real,AbstractArray{#s234,1} where #s234<:Real}",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.translate",
    "category": "method",
    "text": "translate(cluster::AbstractMatrix{<:Real},vector::AbstractVector{<:Real}) -> Matrix{vector|>eltype}\n\nGet the translated cluster of the original one by a vector.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Essentials.Spatials.volume-Tuple{AbstractArray{#s210,1} where #s210<:Real,AbstractArray{#s209,1} where #s209<:Real,AbstractArray{#s208,1} where #s208<:Real}",
    "page": "Spatials",
    "title": "Hamiltonian.Essentials.Spatials.volume",
    "category": "method",
    "text": "volume(v1::AbstractVector{<:Real},v2::AbstractVector{<:Real},v3::AbstractVector{<:Real}) -> Real\n\nGet the volume spanned by three vectors.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Prerequisites.Interfaces.dimension-Tuple{AbstractBond}",
    "page": "Spatials",
    "title": "Hamiltonian.Prerequisites.Interfaces.dimension",
    "category": "method",
    "text": "dimension(bond::AbstractBond) -> Int\ndimension(::Type{<:AbstractBond{R,<:PID,N}}) where {R,N} -> Int\n\nGet the space dimension of a concrete bond.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Prerequisites.Interfaces.dimension-Tuple{AbstractLattice}",
    "page": "Spatials",
    "title": "Hamiltonian.Prerequisites.Interfaces.dimension",
    "category": "method",
    "text": "dimension(lattice::AbstractLattice) -> Int\ndimension(::Type{<:AbstractLattice{<:PID,N}}) where N -> Int\n\nGet the space dimension of the lattice.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Hamiltonian.Prerequisites.Interfaces.rank-Tuple{AbstractBond}",
    "page": "Spatials",
    "title": "Hamiltonian.Prerequisites.Interfaces.rank",
    "category": "method",
    "text": "rank(bond::AbstractBond) -> Int\nrank(::Type{<:AbstractBond{R}}) where R -> Int\n\nGet the rank of a concrete bond.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Base.:==-Tuple{AbstractLattice,AbstractLattice}",
    "page": "Spatials",
    "title": "Base.:==",
    "category": "method",
    "text": "==(lattice1::AbstractLattice,lattice2::AbstractLattice) -> Bool\nisequal(lattice1::AbstractLattice,lattice2::AbstractLattice) -> Bool\n\nOverloaded equivalent operator.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Base.:==-Union{Tuple{R}, Tuple{AbstractBond{R,P,N} where N where P<:PID,AbstractBond{R,P,N} where N where P<:PID}} where R",
    "page": "Spatials",
    "title": "Base.:==",
    "category": "method",
    "text": "==(b1::AbstractBond{R},b1::AbstractBond{R}) where R -> Bool\nisequal(b1::AbstractBond{R},b1::AbstractBond{R}) where R -> Bool\n\nOverloaded equivalent operator.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Base.eltype-Tuple{AbstractBond}",
    "page": "Spatials",
    "title": "Base.eltype",
    "category": "method",
    "text": "eltype(bond::AbstractBond)\neltype(::Type{<:AbstractBond})\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Base.getindex-Tuple{AbstractLattice,RCoordIndex{Int64}}",
    "page": "Spatials",
    "title": "Base.getindex",
    "category": "method",
    "text": "getindex(lattice::AbstractLattice,i::RCoordIndex) -> SVector\ngetindex(lattice::AbstractLattice,i::ICoordIndex) -> SVector\ngetindex(lattice::AbstractLattice,i::PointIndex) -> Point\n\nGet a rcoord, an icoord or a point of a lattice according to the type of the input index.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Base.insert!-Union{Tuple{S}, Tuple{Cylinder,Vararg{S,N} where N}} where S",
    "page": "Spatials",
    "title": "Base.insert!",
    "category": "method",
    "text": "insert!(cylinder::Cylinder,ps::S...;cut::Int=length(cylinder)÷2+1,scopes::Union{<:AbstractVector{S},Nothing}=nothing,coordination::Int=9) where S -> Cylinder\n\nInsert a couple of blocks into a cylinder.\n\nThe position of the cut of the cylinder is specified by the keyword argument cut, which is the center of the cylinder by default. All pids corresponding to a same newly inserted block share the same scope, which is specified by the parameter ps. Optionally, the scopes of the old pids in the cylinder can be replaced if the parameter scopes is assigned other than nothing. Note the length of ps is equal to the number of newly inserted blocks, while that of scopes should be equal to the old length of the cylinder.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Base.iterate",
    "page": "Spatials",
    "title": "Base.iterate",
    "category": "function",
    "text": "iterate(p::Point,state=1)\n\nIterate over the point.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Base.iterate",
    "page": "Spatials",
    "title": "Base.iterate",
    "category": "function",
    "text": "iterate(bond::Bond,state=1)\n\nIterate over the points in a bond.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Base.keytype-Tuple{AbstractLattice}",
    "page": "Spatials",
    "title": "Base.keytype",
    "category": "method",
    "text": "keytype(lattice::AbstractLattice)\nkeytype(::Type{<:AbstractLattice{P}}) where {P<:PID}\n\nGet the pid type of the lattice.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Base.length-Tuple{AbstractBond}",
    "page": "Spatials",
    "title": "Base.length",
    "category": "method",
    "text": "length(bond::AbstractBond) -> Int\nlength(::Type{<:AbstractBond{R}}) where R -> Int\n\nGet the number of points of a bond.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Base.length-Tuple{AbstractLattice}",
    "page": "Spatials",
    "title": "Base.length",
    "category": "method",
    "text": "length(lattice::AbstractLattice) -> Int\n\nGet the number of points contained in a lattice.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Base.reverse-Tuple{Bond}",
    "page": "Spatials",
    "title": "Base.reverse",
    "category": "method",
    "text": "reverse(bond::Bond) -> Bond\n\nGet the reversed bond.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Base.show-Tuple{IO,AbstractLattice}",
    "page": "Spatials",
    "title": "Base.show",
    "category": "method",
    "text": "show(io::IO,lattice::AbstractLattice)\n\nShow a lattice.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Base.show-Tuple{IO,Bond}",
    "page": "Spatials",
    "title": "Base.show",
    "category": "method",
    "text": "show(io::IO,bond::Bond)\n\nShow a bond.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Base.show-Tuple{IO,Point}",
    "page": "Spatials",
    "title": "Base.show",
    "category": "method",
    "text": "show(io::IO,p::Point)\n\nShow a labeled point.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Base.valtype-Tuple{AbstractLattice}",
    "page": "Spatials",
    "title": "Base.valtype",
    "category": "method",
    "text": "valtype(lattice::AbstractLattice)\nvaltype(::Type{<:AbstractLattice{P,N}}) where {P<:PID,N}\n\nGet the point type of the lattice.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Spatials.html#Manul-1",
    "page": "Spatials",
    "title": "Manul",
    "category": "section",
    "text": "Modules=[Spatials]\nOrder=  [:module,:constant,:type,:macro,:function]"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#",
    "page": "Degrees of freedom",
    "title": "Degrees of freedom",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Essentials.DegreesOfFreedom"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Degrees-of-freedom-1",
    "page": "Degrees of freedom",
    "title": "Degrees of freedom",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#IID-and-Internal-1",
    "page": "Degrees of freedom",
    "title": "IID and Internal",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Index-1",
    "page": "Degrees of freedom",
    "title": "Index",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#IDFConfig-1",
    "page": "Degrees of freedom",
    "title": "IDFConfig",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Table-1",
    "page": "Degrees of freedom",
    "title": "Table",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#OID-and-Operator-1",
    "page": "Degrees of freedom",
    "title": "OID and Operator",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Operators-1",
    "page": "Degrees of freedom",
    "title": "Operators",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Hamiltonian.Essentials.DegreesOfFreedom.DirectIndexToTuple",
    "page": "Degrees of freedom",
    "title": "Hamiltonian.Essentials.DegreesOfFreedom.DirectIndexToTuple",
    "category": "type",
    "text": "DirectIndexToTuple\n\nDirect index to tuple.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Hamiltonian.Essentials.DegreesOfFreedom.DirectIndexToTuple-Tuple{Index}",
    "page": "Degrees of freedom",
    "title": "Hamiltonian.Essentials.DegreesOfFreedom.DirectIndexToTuple",
    "category": "method",
    "text": "(indextotuple::DirectIndexToTuple)(index::Index) -> Tuple\n\nConvert an index to tuple directly.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Hamiltonian.Essentials.DegreesOfFreedom.FilteredAttributes",
    "page": "Degrees of freedom",
    "title": "Hamiltonian.Essentials.DegreesOfFreedom.FilteredAttributes",
    "category": "type",
    "text": "FilteredAttributes(::Type{I}) where I<:Index\n\nA method that converts an arbitary index to a tuple, by iterating over the selected attributes in a specific order.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Hamiltonian.Essentials.DegreesOfFreedom.FilteredAttributes-Tuple{Index}",
    "page": "Degrees of freedom",
    "title": "Hamiltonian.Essentials.DegreesOfFreedom.FilteredAttributes",
    "category": "method",
    "text": "(indextotuple::FilteredAttributes)(index::Index) -> Tuple\n\nConvert an index to tuple by the \"filtered attributes\" method.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Hamiltonian.Essentials.DegreesOfFreedom.IDFConfig",
    "page": "Degrees of freedom",
    "title": "Hamiltonian.Essentials.DegreesOfFreedom.IDFConfig",
    "category": "type",
    "text": "IDFConfig{I}(map::Function,pids::Union{AbstractVector{<:PID},Tuple{}}=()) where I<:Internal\n\nConfiguration of the internal degrees of freedom at a lattice.\n\nmap maps a PID to an Internal.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Hamiltonian.Essentials.DegreesOfFreedom.IID",
    "page": "Degrees of freedom",
    "title": "Hamiltonian.Essentials.DegreesOfFreedom.IID",
    "category": "type",
    "text": "IID\n\nThe id of an internal degree of freedom.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Hamiltonian.Essentials.DegreesOfFreedom.Index",
    "page": "Degrees of freedom",
    "title": "Hamiltonian.Essentials.DegreesOfFreedom.Index",
    "category": "type",
    "text": "Index{P,I}\n\nThe complete index of a degree of freedom, which consist of the spatial part and the internal part.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Hamiltonian.Essentials.DegreesOfFreedom.IndexToTuple",
    "page": "Degrees of freedom",
    "title": "Hamiltonian.Essentials.DegreesOfFreedom.IndexToTuple",
    "category": "type",
    "text": "IndexToTuple\n\nThe rules for converting an index to a tuple.\n\nAs a function, every instance should accept only one positional argument, i.e. the index to be converted to a tuple.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Hamiltonian.Essentials.DegreesOfFreedom.Internal",
    "page": "Degrees of freedom",
    "title": "Hamiltonian.Essentials.DegreesOfFreedom.Internal",
    "category": "type",
    "text": "Internal\n\nThe whole internal degrees of freedom at a single point.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Hamiltonian.Essentials.DegreesOfFreedom.OID",
    "page": "Degrees of freedom",
    "title": "Hamiltonian.Essentials.DegreesOfFreedom.OID",
    "category": "type",
    "text": "OID(index::Index,::Nothing,::Nothing,seq::Union{Nothing,Int})\nOID(index::Index,rcoord::SVector{N,Float},icoord::SVector{N,Float},seq::Union{Nothing,Int}) where N\nOID(index::Index,rcoord::Vector{Float},icoord::Vector{Float},seq::Union{Nothing,Int})\nOID(index::Index;rcoord::Union{Nothing,SVector,Vector{Float}}=nothing,icoord::Union{Nothing,SVector,Vector{Float}}=nothing,seq::Union{Nothing,Int}=nothing)\n\nOperator id.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Hamiltonian.Essentials.DegreesOfFreedom.Operator",
    "page": "Degrees of freedom",
    "title": "Hamiltonian.Essentials.DegreesOfFreedom.Operator",
    "category": "type",
    "text": "Operator{N,V<:Number,I<:ID{<:NTuple{N,OID}}} <: Element{N,V,I}\n\nAbstract type for an operator.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Hamiltonian.Essentials.DegreesOfFreedom.Operators",
    "page": "Degrees of freedom",
    "title": "Hamiltonian.Essentials.DegreesOfFreedom.Operators",
    "category": "type",
    "text": "Operators(opts::Operator...)\n\nA set of operators.\n\nType alias of Operators{I<:ID,O<:Operator}=Elements{I,O}.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Hamiltonian.Essentials.DegreesOfFreedom.Table",
    "page": "Degrees of freedom",
    "title": "Hamiltonian.Essentials.DegreesOfFreedom.Table",
    "category": "type",
    "text": "Table{I<:Index} <: AbstractDict{I,Int}\n\nIndex-sequence table. Alias for Dict{I<:Index,Int}.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Hamiltonian.Essentials.DegreesOfFreedom.Table-Tuple{AbstractArray{#s233,1} where #s233<:Index}",
    "page": "Degrees of freedom",
    "title": "Hamiltonian.Essentials.DegreesOfFreedom.Table",
    "category": "method",
    "text": "Table(indices::AbstractVector{<:Index};by::IndexToTuple=directindextotuple) -> Table\n\nConvert an sequence of indices to the corresponding index-sequence table.\n\nThe input indices will be converted to tuples by the by function with the duplicates removed. The resulting unique tuples are sorted, which determines the sequence of the input indices. Note that two indices have the same sequence if their converted tupels are equal to each other.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Hamiltonian.Essentials.DegreesOfFreedom.Table-Tuple{IDFConfig}",
    "page": "Degrees of freedom",
    "title": "Hamiltonian.Essentials.DegreesOfFreedom.Table",
    "category": "method",
    "text": "Table(config::IDFConfig;by::IndexToTuple=directindextotuple) -> Table\n\nGet the index-sequence table of the whole internal Hilbert spaces at a lattice.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Hamiltonian.Essentials.DegreesOfFreedom.directindextotuple",
    "page": "Degrees of freedom",
    "title": "Hamiltonian.Essentials.DegreesOfFreedom.directindextotuple",
    "category": "function",
    "text": "directindextotuple\n\nIndicate that the conversion from an index to a tuple is direct.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Hamiltonian.Essentials.DegreesOfFreedom.iid-Tuple{Index}",
    "page": "Degrees of freedom",
    "title": "Hamiltonian.Essentials.DegreesOfFreedom.iid",
    "category": "method",
    "text": "iid(index::Index) -> IID\n\nGet the internal part of an index.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Hamiltonian.Essentials.DegreesOfFreedom.iidtype-Tuple{Index}",
    "page": "Degrees of freedom",
    "title": "Hamiltonian.Essentials.DegreesOfFreedom.iidtype",
    "category": "method",
    "text": "iidtype(index::Index)\niidtype(::Type{<:Index{P,I}}) where {P,I}\n\nGet the type of the internal part of an index.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Hamiltonian.Essentials.DegreesOfFreedom.isHermitian-Tuple{Dict{I,O} where O<:Operator where I<:Hamiltonian.Mathematics.AlgebraOverFields.ID}",
    "page": "Degrees of freedom",
    "title": "Hamiltonian.Essentials.DegreesOfFreedom.isHermitian",
    "category": "method",
    "text": "isHermitian(opts::Operators) -> Bool\n\nJudge whether a set of operators as a whole is Hermitian.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Hamiltonian.Essentials.DegreesOfFreedom.isHermitian-Tuple{Operator}",
    "page": "Degrees of freedom",
    "title": "Hamiltonian.Essentials.DegreesOfFreedom.isHermitian",
    "category": "method",
    "text": "isHermitian(opt::Operator) -> Bool\n\nJudge whether an operator is Hermitian.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Hamiltonian.Essentials.DegreesOfFreedom.isHermitian-Union{Tuple{ID{#s236} where #s236<:Tuple{Vararg{OID,N}}}, Tuple{N}} where N",
    "page": "Degrees of freedom",
    "title": "Hamiltonian.Essentials.DegreesOfFreedom.isHermitian",
    "category": "method",
    "text": "isHermitian(oid::ID{<:NTuple{N,OID}}) where N -> Bool\n\nJudge whether an operator id is Hermitian.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Hamiltonian.Essentials.DegreesOfFreedom.oidtype",
    "page": "Degrees of freedom",
    "title": "Hamiltonian.Essentials.DegreesOfFreedom.oidtype",
    "category": "function",
    "text": "oidtype\n\nGet the compatible oid type.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Hamiltonian.Essentials.DegreesOfFreedom.otype",
    "page": "Degrees of freedom",
    "title": "Hamiltonian.Essentials.DegreesOfFreedom.otype",
    "category": "function",
    "text": "otype\n\nGet the compatible operator type from a term type, a bond type and a table type.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Hamiltonian.Essentials.DegreesOfFreedom.pid-Tuple{Index}",
    "page": "Degrees of freedom",
    "title": "Hamiltonian.Essentials.DegreesOfFreedom.pid",
    "category": "method",
    "text": "pid(index::Index) -> PID\n\nGet the spatial part of an index.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Hamiltonian.Essentials.Spatials.icoord-Tuple{Operator{1,V,I} where I<:(Hamiltonian.Mathematics.AlgebraOverFields.ID{#s254} where #s254<:Tuple{OID}) where V<:Number}",
    "page": "Degrees of freedom",
    "title": "Hamiltonian.Essentials.Spatials.icoord",
    "category": "method",
    "text": "icoord(opt::Operator{1}) -> SVector\nicoord(opt::Operator{2}) -> SVector\n\nGet the whole icoord of an operator.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Hamiltonian.Essentials.Spatials.pidtype-Tuple{Index}",
    "page": "Degrees of freedom",
    "title": "Hamiltonian.Essentials.Spatials.pidtype",
    "category": "method",
    "text": "pidtype(index::Index)\npidtype(::Type{<:Index{P,I}}) where {P,I}\n\nGet the type of the spatial part of an index.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Hamiltonian.Essentials.Spatials.rcoord-Tuple{Operator{1,V,I} where I<:(Hamiltonian.Mathematics.AlgebraOverFields.ID{#s254} where #s254<:Tuple{OID}) where V<:Number}",
    "page": "Degrees of freedom",
    "title": "Hamiltonian.Essentials.Spatials.rcoord",
    "category": "method",
    "text": "rcoord(opt::Operator{1}) -> SVector\nrcoord(opt::Operator{2}) -> SVector\n\nGet the whole rcoord of an operator.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Core.Type-Tuple{PID,IID}",
    "page": "Degrees of freedom",
    "title": "Core.Type",
    "category": "method",
    "text": "(INDEX::Type{<:Index})(pid::PID,iid::IID) -> INDEX\n\nGet the corresponding index from a pid and an iid.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Base.adjoint-Tuple{Dict{I,O} where O<:Operator where I<:Hamiltonian.Mathematics.AlgebraOverFields.ID}",
    "page": "Degrees of freedom",
    "title": "Base.adjoint",
    "category": "method",
    "text": "adjoint(opts::Operators) -> Operators\n\nGet the adjoint of a set of operators.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Base.adjoint-Tuple{Index}",
    "page": "Degrees of freedom",
    "title": "Base.adjoint",
    "category": "method",
    "text": "adjoint(index::Index) -> typeof(index)\n\nGet the adjoint of an index.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Base.adjoint-Tuple{OID}",
    "page": "Degrees of freedom",
    "title": "Base.adjoint",
    "category": "method",
    "text": "adjoint(oid::OID) -> typeof(oid)\nadjoint(oid::ID{<:NTuple{N,OID}}) where N -> typeof(oid)\n\nGet the adjoint of an operator id.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Base.adjoint-Tuple{Operator}",
    "page": "Degrees of freedom",
    "title": "Base.adjoint",
    "category": "method",
    "text": "adjoint(opt::Operator{N}) where N -> Operator\n\nGet the adjoint of an operator.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Base.filter-Tuple{Function,FilteredAttributes}",
    "page": "Degrees of freedom",
    "title": "Base.filter",
    "category": "method",
    "text": "filter(f::Function,indextotuple::FilteredAttributes)\n\nFilter the attributes of a \"filtered attributes\" method.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Base.length-Tuple{FilteredAttributes}",
    "page": "Degrees of freedom",
    "title": "Base.length",
    "category": "method",
    "text": "length(indextotuple::FilteredAttributes) -> Int\nlength(::Type{<:FilteredAttributes{N}}) where N -> Int\n\nGet the length of the filtered attributes.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Base.replace!-Tuple{IDFConfig,Vararg{PID,N} where N}",
    "page": "Degrees of freedom",
    "title": "Base.replace!",
    "category": "method",
    "text": "replace!(config::IDFConfig,pids::PID...) -> IDFConfig\n\nReset the idfconfig with new pids.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Base.reverse-Tuple{Dict{I,Int64} where I<:Index}",
    "page": "Degrees of freedom",
    "title": "Base.reverse",
    "category": "method",
    "text": "reverse(table::Table) -> Dict{Int,Set{<:Index}}\n\nConvert an index-sequence table to a sequence-indices table.\n\nSince different indices may correspond to the same sequence, the reverse is a one-to-many map.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Base.show-Tuple{IO,Internal}",
    "page": "Degrees of freedom",
    "title": "Base.show",
    "category": "method",
    "text": "show(io::IO,i::Internal)\n\nShow an internal.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Base.show-Tuple{IO,OID}",
    "page": "Degrees of freedom",
    "title": "Base.show",
    "category": "method",
    "text": "show(io::IO,oid::OID)\n\nShow an operator id.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Base.show-Tuple{IO,Operator}",
    "page": "Degrees of freedom",
    "title": "Base.show",
    "category": "method",
    "text": "show(io::IO,opt::Operator)\n\nShow an operator.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Base.union-Tuple{Vararg{Dict{I,Int64} where I<:Index,N} where N}",
    "page": "Degrees of freedom",
    "title": "Base.union",
    "category": "method",
    "text": "union(tables::Table...;by::IndexToTuple=directindextotuple) -> Table\n\nUnite several index-sequence tables.\n\nSee Table for more details.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Base.union-Union{Tuple{I}, Tuple{P}, Tuple{Type{P},Type{I}}} where I<:IID where P<:PID",
    "page": "Degrees of freedom",
    "title": "Base.union",
    "category": "method",
    "text": "union(::Type{P},::Type{I}) where {P<:PID,I<:IID}\n\nCombine a concrete PID type and a concrete IID type to a concrete Index type.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/DegreesOfFreedom.html#Manul-1",
    "page": "Degrees of freedom",
    "title": "Manul",
    "category": "section",
    "text": "Modules=[DegreesOfFreedom]\nOrder=  [:module,:constant,:type,:macro,:function]"
},

{
    "location": "man/Essentials/Terms.html#",
    "page": "Terms",
    "title": "Terms",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Essentials.Termspush!(LOAD_PATH,\"../../../../src/\")\nusing Hamiltonian"
},

{
    "location": "man/Essentials/Terms.html#Terms-1",
    "page": "Terms",
    "title": "Terms",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/Terms.html#Subscript-and-Subscripts-1",
    "page": "Terms",
    "title": "Subscript and Subscripts",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/Terms.html#Coupling-and-Couplings-1",
    "page": "Terms",
    "title": "Coupling and Couplings",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/Terms.html#Term-1",
    "page": "Terms",
    "title": "Term",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/Terms.html#Hamiltonian.Essentials.Terms.Coupling",
    "page": "Terms",
    "title": "Hamiltonian.Essentials.Terms.Coupling",
    "category": "type",
    "text": "Coupling{N,V<:Number,I<:ID{<:NTuple{N,SimpleID}}} <: Element{N,V,I}\n\nThe coupling intra/inter interanl degrees of freedom at different lattice points.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Hamiltonian.Essentials.Terms.Couplings",
    "page": "Terms",
    "title": "Hamiltonian.Essentials.Terms.Couplings",
    "category": "type",
    "text": "Couplings{I<:ID,C<:Coupling} <: AbstractDict{I,C}\n\nA pack of couplings intra/inter interanl degrees of freedom at different lattice points.\n\nAlias for Elements{I,C}.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Hamiltonian.Essentials.Terms.Subscript",
    "page": "Terms",
    "title": "Hamiltonian.Essentials.Terms.Subscript",
    "category": "type",
    "text": "Subscript(  ipattern::NTuple{N1,Any},\n            opattern::NTuple{N2,Any},\n            mapping::Union{Function,Nothing}=nothing,\n            constrain::Union{Function,Nothing}=nothing,\n            identifier::Union{Symbol,Char}=wildcard\n            ) where {N1,N2}\nSubscript{N}() where N\nSubscript(opattern::NTuple{N,Int}) where N\n\nThe subscripts of some orbital/spin degrees of freedom.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Hamiltonian.Essentials.Terms.Subscript-Tuple{}",
    "page": "Terms",
    "title": "Hamiltonian.Essentials.Terms.Subscript",
    "category": "method",
    "text": "(subscript::Subscript{N})(::Val{\'M\'},values::Vararg{Int,N}) where N -> NTuple{dimension(subscript),Int}\n(subscript::Subscript{N})(::Val{\'C\'},values::Vararg{Int,N}) where N -> Bool\n\nConstruct the subscripts from a set of independent variables.\nJudge whether a set of independent variables are valid to construct the subscripts.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Hamiltonian.Essentials.Terms.Subscripts",
    "page": "Terms",
    "title": "Hamiltonian.Essentials.Terms.Subscripts",
    "category": "type",
    "text": "Subscripts(contents::Subscript...)\n\nA complete set of all the independent subscripts of the orbital/spin degrees of freedom.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Hamiltonian.Essentials.Terms.Subscripts-Tuple{}",
    "page": "Terms",
    "title": "Hamiltonian.Essentials.Terms.Subscripts",
    "category": "method",
    "text": "(subscripts::Subscripts)(::Val{\'M\'},values::NTuple{N,Int}) where N -> NTuple{dimension(subscripts),Int}\n(subscripts::Subscripts)(::Val{\'C\'},values::NTuple{N,Int}) where N -> Bool\n\nConstruct the complete set of subscripts from a complete set of independent variables.\nJudge whether a complete set of independent variables are valid to construct the complete subscripts.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Hamiltonian.Essentials.Terms.Term",
    "page": "Terms",
    "title": "Hamiltonian.Essentials.Terms.Term",
    "category": "type",
    "text": "Term{ST,SP}(id::Symbol,value::Number,neighbor::Any,couplings::TermCouplings,amplitude::TermAmplitude,modulate::Union{TermModulate,Nothing},factor::Number) where {ST,SP}\nTerm{ST,SP}(id::Symbol,value::Number,neighbor::Any;\n            couplings::Union{Tuple{<:Tuple{Vararg{Couplings}},<:Function},Coupling,Couplings,TermCouplings},\n            amplitude::Union{Function,Nothing}=nothing,\n            modulate::Union{Function,Bool}=false,\n            factor::Number=1\n            ) where {ST,SP}\n\nA term of a quantum lattice system.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Hamiltonian.Essentials.Terms.TermAmplitude",
    "page": "Terms",
    "title": "Hamiltonian.Essentials.Terms.TermAmplitude",
    "category": "type",
    "text": "TermAmplitude(amplitude::Union{Function,Nothing}=nothing)\n\nThe function for the amplitude of a term.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Hamiltonian.Essentials.Terms.TermCouplings",
    "page": "Terms",
    "title": "Hamiltonian.Essentials.Terms.TermCouplings",
    "category": "type",
    "text": "TermCouplings(candidate::Coupling)\nTermCouplings(candidate::Couplings)\nTermCouplings(contents::Tuple{<:Tuple{Vararg{Couplings}},<:Function})\nTermCouplings(candidates::NTuple{N,<:Couplings},choice::Function) where N\n\nThe function for the couplings of a term.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Hamiltonian.Essentials.Terms.TermFunction",
    "page": "Terms",
    "title": "Hamiltonian.Essentials.Terms.TermFunction",
    "category": "type",
    "text": "TermFunction <: Function\n\nAbstract type for concrete term functions.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Hamiltonian.Essentials.Terms.TermModulate",
    "page": "Terms",
    "title": "Hamiltonian.Essentials.Terms.TermModulate",
    "category": "type",
    "text": "TermModulate(id::Symbol,modulate::Union{Function,Nothing}=nothing)\nTermModulate(id::Symbol,modulate::Bool)\n\nThe function for the modulation of a term.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Hamiltonian.Essentials.Terms.@subscript",
    "page": "Terms",
    "title": "Hamiltonian.Essentials.Terms.@subscript",
    "category": "macro",
    "text": "@subscript expr::Expr with constrain::Expr -> Subscript\n\nConstruct a subscript from a map and optionally with a constrain.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Hamiltonian.Essentials.Terms.abbr-Tuple{Term}",
    "page": "Terms",
    "title": "Hamiltonian.Essentials.Terms.abbr",
    "category": "method",
    "text": "abbr(term::Term) -> Symbol\nabbr(::Type{<:Term}) -> Symbol\n\nGet the abbreviation of the species of a term.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Hamiltonian.Essentials.Terms.species-Tuple{Term}",
    "page": "Terms",
    "title": "Hamiltonian.Essentials.Terms.species",
    "category": "method",
    "text": "species(term::Term) -> Symbol\nspecies(::Type{<:Term{ST,SP}}) where {ST,SP} -> Symbol\n\nGet the species of a term.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Hamiltonian.Essentials.Terms.statistics-Tuple{Term}",
    "page": "Terms",
    "title": "Hamiltonian.Essentials.Terms.statistics",
    "category": "method",
    "text": "statistics(term::Term) -> Char\nstatistics(::Type{<:Term{ST}}) where ST -> Char\n\nGet the statistics of a term.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Hamiltonian.Prerequisites.Interfaces.expand",
    "page": "Terms",
    "title": "Hamiltonian.Prerequisites.Interfaces.expand",
    "category": "function",
    "text": "expand(otype::Type{<:Operator},term::Term,bond::AbstractBond,config::IDFConfig,table::Union{Table,Nothing}=nothing,half::Union{Bool,Nothing}=nothing)\n\nExpand the operators of a term on a bond with a given config.\n\nThe half parameter determines the behavior of generating operators, which falls into the following three categories\n\nfalse: no extra operations on the generated operators\ntrue: an extra multiplication by 0.5 with the generated operators\nnothing: \"Hermitian half\" of the generated operators\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Hamiltonian.Prerequisites.Interfaces.expand-Union{Tuple{N}, Tuple{Subscripts,Tuple{Vararg{Int64,N}}}} where N",
    "page": "Terms",
    "title": "Hamiltonian.Prerequisites.Interfaces.expand",
    "category": "method",
    "text": "expand(subscripts::Subscripts,dimensions::NTuple{N,Int}) where N -> SbExpand\n\nExpand a complete set of subscripts with a given set of variable ranges.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Hamiltonian.Prerequisites.Interfaces.rank-Tuple{Subscripts,Int64}",
    "page": "Terms",
    "title": "Hamiltonian.Prerequisites.Interfaces.rank",
    "category": "method",
    "text": "rank(subscripts::Subscripts,i::Int) -> Int\nrank(::Type{<:Subscripts{T}},i::Int) where T -> Int\n\nGet the number of the independent variables of a component of the complete subscript set.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Hamiltonian.Prerequisites.Interfaces.rank-Tuple{Subscripts}",
    "page": "Terms",
    "title": "Hamiltonian.Prerequisites.Interfaces.rank",
    "category": "method",
    "text": "rank(subscripts::Subscripts) -> Int\nrank(::Type{S}) where S<:Subscripts -> Int\n\nGet the total number of the independent variables of the complete subscript set.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Hamiltonian.Prerequisites.Interfaces.rank-Tuple{Subscript}",
    "page": "Terms",
    "title": "Hamiltonian.Prerequisites.Interfaces.rank",
    "category": "method",
    "text": "rank(subscript::Subscript) -> Int\nrank(::Type{<:Subscript{N}}) where N -> Int\n\nGet the number of the independent variables that are used to describe the subscripts of some orbital/spin degrees of freedom.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Hamiltonian.Prerequisites.Interfaces.rank-Tuple{TermCouplings}",
    "page": "Terms",
    "title": "Hamiltonian.Prerequisites.Interfaces.rank",
    "category": "method",
    "text": "rank(tcs::TermCouplings) -> Int\nrank(TCS::Type{<:TermCouplings}) -> Int\n\nGet the rank of the couplings it contained.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Hamiltonian.Prerequisites.Interfaces.rank-Tuple{Term}",
    "page": "Terms",
    "title": "Hamiltonian.Prerequisites.Interfaces.rank",
    "category": "method",
    "text": "rank(term::Term) -> Int\nrank(::Type{T}) where T<:Term -> Int\n\nGet the rank of a term.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Hamiltonian.Prerequisites.Interfaces.update!-Tuple{Term,Vararg{Any,N} where N}",
    "page": "Terms",
    "title": "Hamiltonian.Prerequisites.Interfaces.update!",
    "category": "method",
    "text": "update!(term::Term,args...;kwargs...) -> Term\n\nUpdate the value of a term by its modulate function.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Base.:*-Tuple{Subscript,Subscript}",
    "page": "Terms",
    "title": "Base.:*",
    "category": "method",
    "text": "*(sub1::Subscript,sub2::Subscript) -> Subscripts\n*(subs::Subscripts,sub::Subscript) -> Subscripts\n*(sub::Subscript,subs::Subscripts) -> Subscripts\n*(subs1::Subscripts,subs2::Subscripts) -> Subscripts\n\nGet the multiplication between subscripts or complete sets of subscripts.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Base.:+-Tuple{Term}",
    "page": "Terms",
    "title": "Base.:+",
    "category": "method",
    "text": "+(term::Term) -> Term\n-(term::Term) -> Term\n*(term::Term,factor::Number) -> Term\n*(factor::Number,term::Term) -> Term\n/(term::Term,factor::Number) -> Term\n\nAllowed arithmetic operations for a term.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Base.:==-Tuple{Subscript,Subscript}",
    "page": "Terms",
    "title": "Base.:==",
    "category": "method",
    "text": "==(sub1::Subscript,sub2::Subscript) -> Bool\nisequal(sub1::Subscript,sub2::Subscript) -> Bool\n\nJudge whether two subscripts are equivalent to each other.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Base.:==-Tuple{Term,Term}",
    "page": "Terms",
    "title": "Base.:==",
    "category": "method",
    "text": "==(term1::Term,term2::Term) -> Bool\nisequal(term1::Term,term2::Term) -> Bool\n\nJudge whether two terms are equivalent to each other.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Base.:==-Tuple{TermFunction,TermFunction}",
    "page": "Terms",
    "title": "Base.:==",
    "category": "method",
    "text": "==(tf1::TermFunction,tf2::TermFunction) -> Bool\nisequal(tf1::TermFunction,tf2::TermFunction) -> Bool\n\nJudge whether two concrete term functions are equivalent to each other.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Base.one-Tuple{Term}",
    "page": "Terms",
    "title": "Base.one",
    "category": "method",
    "text": "one(term::Term) -> Term\n\nGet a unit term.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Base.replace-Union{Tuple{Term{ST,SP,V,N,C,A,M} where M<:Union{Nothing, TermModulate} where A<:TermAmplitude where C<:TermCouplings where N where V<:Number}, Tuple{SP}, Tuple{ST}} where SP where ST",
    "page": "Terms",
    "title": "Base.replace",
    "category": "method",
    "text": "replace(term::Term{ST,SP};kwargs...) where {ST,SP} -> Term\n\nReplace some attributes of a term with key word arguments.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Base.repr-Tuple{Term,AbstractBond,IDFConfig}",
    "page": "Terms",
    "title": "Base.repr",
    "category": "method",
    "text": "repr(term::Term,bond::AbstractBond,config::IDFConfig) -> String\n\nGet the repr representation of a term on a bond with a given config.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Base.show-Tuple{IO,Subscript}",
    "page": "Terms",
    "title": "Base.show",
    "category": "method",
    "text": "show(io::IO,subscript::Subscript)\n\nShow a subscript.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Base.show-Tuple{IO,Term}",
    "page": "Terms",
    "title": "Base.show",
    "category": "method",
    "text": "show(io::IO,term::Term)\n\nShow a term.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Base.valtype-Tuple{Term}",
    "page": "Terms",
    "title": "Base.valtype",
    "category": "method",
    "text": "valtype(term::Term)\nvaltype(::Type{<:Term{ST,SP,V}}) where {ST,SP,V<:Number}\n\nGet the value type of a term.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Base.zero-Tuple{Term}",
    "page": "Terms",
    "title": "Base.zero",
    "category": "method",
    "text": "zero(term::Term) -> Term\n\nGet a zero term.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Hamiltonian.Prerequisites.Interfaces.dimension-Tuple{Subscripts,Int64}",
    "page": "Terms",
    "title": "Hamiltonian.Prerequisites.Interfaces.dimension",
    "category": "method",
    "text": "dimension(subscripts::Subscripts,i::Int) -> Int\ndimension(::Type{<:Subscripts{T}},i::Int) where T -> Int\n\nGet the total number of the whole variables of a component of the complete subscript set.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Hamiltonian.Prerequisites.Interfaces.dimension-Tuple{Subscripts}",
    "page": "Terms",
    "title": "Hamiltonian.Prerequisites.Interfaces.dimension",
    "category": "method",
    "text": "dimension(subscripts::Subscripts) -> Int\ndimension(::Type{S}) where S<:Subscripts -> Int\n\nGet the total number of the whole variables of the complete subscript set.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Hamiltonian.Prerequisites.Interfaces.dimension-Tuple{Subscript}",
    "page": "Terms",
    "title": "Hamiltonian.Prerequisites.Interfaces.dimension",
    "category": "method",
    "text": "dimension(subscript::Subscript) -> Int\ndimension(::Type{<:Subscript{N1,N2}}) where {N1,N2} -> Int\n\nGet the number of the whole variables that are used to describe the subscripts of some orbital/spin degrees of freedom.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/Terms.html#Manul-1",
    "page": "Terms",
    "title": "Manul",
    "category": "section",
    "text": "Modules=[Terms]\nOrder=  [:module,:constant,:type,:macro,:function]"
},

{
    "location": "man/Essentials/FockPackage.html#",
    "page": "Fock package",
    "title": "Fock package",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Essentials.FockPackage"
},

{
    "location": "man/Essentials/FockPackage.html#Fock-package-1",
    "page": "Fock package",
    "title": "Fock package",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/FockPackage.html#Fock-degrees-of-freedom-1",
    "page": "Fock package",
    "title": "Fock degrees of freedom",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/FockPackage.html#FID-and-Fock-1",
    "page": "Fock package",
    "title": "FID and Fock",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/FockPackage.html#FIndex-1",
    "page": "Fock package",
    "title": "FIndex",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/FockPackage.html#Predefined-Fock-operators-1",
    "page": "Fock package",
    "title": "Predefined Fock operators",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/FockPackage.html#Fock-terms-1",
    "page": "Fock package",
    "title": "Fock terms",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/FockPackage.html#FCID-and-FockCoupling-1",
    "page": "Fock package",
    "title": "FCID and FockCoupling",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/FockPackage.html#Predefined-Fock-couplings-1",
    "page": "Fock package",
    "title": "Predefined Fock couplings",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/FockPackage.html#Predefined-Fock-terms-1",
    "page": "Fock package",
    "title": "Predefined Fock terms",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.ANNIHILATION",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.ANNIHILATION",
    "category": "constant",
    "text": "ANNIHILATION\n\nIndicate that the nambu index is ANNIHILATION.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.CREATION",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.CREATION",
    "category": "constant",
    "text": "CREATION\n\nIndicate that the nambu index is CREATION.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.BOperator",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.BOperator",
    "category": "type",
    "text": "BOperator(value::Number,id::ID{<:NTuple{N,OID}}) where N\n\nBosonic Fock operator.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.Coulomb",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.Coulomb",
    "category": "type",
    "text": "Coulomb{ST}(    id::Symbol,value::Number;\n                neighbor::Int=1,\n                couplings::Union{Tuple{<:Tuple{Vararg{Couplings}},<:Function},Coupling,Couplings,Nothing}=nothing,\n                amplitude::Union{Function,Nothing}=nothing,\n                modulate::Union{Function,Bool}=false,\n                factor::Number=1\n                ) where {ST}\n\nCoulomb term.\n\nType alias for Term{Statistics,:Coulomb,<:Number,Int,<:TermCouplings,<:TermAmplitude,<:Union{TermModulate,Nothing}}.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.FCID",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.FCID",
    "category": "type",
    "text": "FCID(;center=wildcard,atom=wildcard,orbital=wildcard,spin=wildcard,nambu=wildcard,obsub=wildcard,spsub=wildcard)\n\nThe id of a Fock coupling.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.FID",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.FID",
    "category": "type",
    "text": "FID <: IID\n\nThe Fock id.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.FID-Tuple{}",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.FID",
    "category": "method",
    "text": "FID(;orbital::Int=1,spin::Int=1,nambu::Int=ANNIHILATION)\n\nCreate a Fock id.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.FIndex",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.FIndex",
    "category": "type",
    "text": "FIndex{S} <: Index{PID{S},FID}\n\nThe Fock index.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.FOperator",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.FOperator",
    "category": "type",
    "text": "FOperator(value::Number,id::ID{<:NTuple{N,OID}}) where N\n\nFermionic Fock operator.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.Fock",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.Fock",
    "category": "type",
    "text": "Fock <: Internal{FID}\n\nThe Fock interanl degrees of freedom.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.Fock-Tuple{}",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.Fock",
    "category": "method",
    "text": "Fock(;atom::Int=1,norbital::Int=1,nspin::Int=2,nnambu::Int=2)\n\nConstruct a Fock degrees of freedom.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.FockCoupling",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.FockCoupling",
    "category": "type",
    "text": "FockCoupling(value::Number,id::ID{<:NTuple{N,FCID}},obsubscripts::Subscripts,spsubscripts::Subscripts) where N\nFockCoupling{N}(value::Number=1;\n                centers::Union{NTuple{N,Int},Nothing}=nothing,\n                atoms::Union{NTuple{N,Int},Nothing}=nothing,\n                orbitals::Union{NTuple{N,Int},Subscript,Nothing}=nothing,\n                spins::Union{NTuple{N,Int},Subscript,Nothing}=nothing,\n                nambus::Union{NTuple{N,Int},Nothing}=nothing) where N\n\nFock coupling.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.FockOperator",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.FockOperator",
    "category": "type",
    "text": "FockOperator{N,V<:Number,I<:ID{<:NTuple{N,OID}}} <: Operator{N,V,I}\n\nAbstract type for all Fock operators.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.Hopping",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.Hopping",
    "category": "type",
    "text": "Hopping{ST}(id::Symbol,value::Number;\n            neighbor::Int=1,\n            couplings::Union{Tuple{<:Tuple{Vararg{Couplings}},<:Function},Coupling,Couplings,Nothing}=nothing,\n            amplitude::Union{Function,Nothing}=nothing,\n            modulate::Union{Function,Bool}=false,\n            factor::Number=1\n            ) where {ST}\n\nHopping term.\n\nType alias for Term{Statistics,:Hopping,<:Number,Int,<:TermCouplings,<:TermAmplitude,<:Union{TermModulate,Nothing}}.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.Hubbard",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.Hubbard",
    "category": "type",
    "text": "Hubbard{ST}(id::Symbol,value::Real;\n            amplitude::Union{Function,Nothing}=nothing,\n            modulate::Union{Function,Bool}=false,\n            factor::Number=1\n            ) where {ST}\n\nHubbard term.\n\nType alias for Term{Statistics,:Hubbard,<:Number,Int,<:TermCouplings,<:TermAmplitude,<:Union{TermModulate,Nothing}}.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.InterOrbitalInterSpin",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.InterOrbitalInterSpin",
    "category": "type",
    "text": "InterOrbitalInterSpin{ST}(  id::Symbol,value::Real;\n                            amplitude::Union{Function,Nothing}=nothing,\n                            modulate::Union{Function,Bool}=false,\n                            factor::Number=1\n                            ) where {ST}\n\nInterorbital-interspin term.\n\nType alias for Term{Statistics,:InterOrbitalInterSpin,<:Number,Int,<:TermCouplings,<:TermAmplitude,<:Union{TermModulate,Nothing}}.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.InterOrbitalIntraSpin",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.InterOrbitalIntraSpin",
    "category": "type",
    "text": "InterOrbitalIntraSpin{ST}(  id::Symbol,value::Real;\n                            amplitude::Union{Function,Nothing}=nothing,\n                            modulate::Union{Function,Bool}=false,\n                            factor::Number=1\n                            ) where {ST}\n\nInterorbital-intraspin term.\n\nType alias for Term{Statistics,:InterOrbitalIntraSpin,<:Number,Int,<:TermCouplings,<:TermAmplitude,<:Union{TermModulate,Nothing}}.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.Onsite",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.Onsite",
    "category": "type",
    "text": "Onsite{ST}( id::Symbol,value::Number;\n            couplings::Union{Tuple{<:Tuple{Vararg{Couplings}},<:Function},Coupling,Couplings,Nothing}=nothing,\n            amplitude::Union{Function,Nothing}=nothing,\n            modulate::Union{Function,Bool}=false,\n            factor::Number=1\n            ) where {ST}\n\nOnsite term.\n\nType alias for Term{Statistics,:Onsite,<:Number,Int,<:TermCouplings,<:TermAmplitude,<:Union{TermModulate,Nothing}}.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.PairHopping",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.PairHopping",
    "category": "type",
    "text": "PairHopping{ST}(    id::Symbol,value::Real;\n                    amplitude::Union{Function,Nothing}=nothing,\n                    modulate::Union{Function,Bool}=false,\n                    factor::Number=1\n                    ) where {ST}\n\nPair-hopping term.\n\nType alias for Term{Statistics,:PairHopping,<:Number,Int,<:TermCouplings,<:TermAmplitude,<:Union{TermModulate,Nothing}}.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.Pairing",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.Pairing",
    "category": "type",
    "text": "Pairing{ST}(id::Symbol,value::Number;\n            neighbor::Int=0,\n            couplings::Union{Tuple{<:Tuple{Vararg{Couplings}},<:Function},Coupling,Couplings,Nothing}=nothing,\n            amplitude::Union{Function,Nothing}=nothing,\n            modulate::Union{Function,Bool}=false,\n            factor::Number=1\n            ) where {ST}\n\nPairing term.\n\nType alias for Term{Statistics,:Pairing,<:Number,Int,<:TermCouplings,<:TermAmplitude,<:Union{TermModulate,Nothing}}.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.SpinFlip",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.SpinFlip",
    "category": "type",
    "text": "SpinFlip{ST}(   id::Symbol,value::Real;\n                amplitude::Union{Function,Nothing}=nothing,\n                modulate::Union{Function,Bool}=false,\n                factor::Number=1\n                ) where {ST}\n\nSpin-flip term.\n\nType alias for Term{Statistics,:SpinFlip,<:Number,Int,<:TermCouplings,<:TermAmplitude,<:Union{TermModulate,Nothing}}.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.DegreesOfFreedom.oidtype-Tuple{Val{:Fock},Type{#s237} where #s237<:AbstractBond,Type{Nothing}}",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.DegreesOfFreedom.oidtype",
    "category": "method",
    "text": "oidtype(::Val{:Fock},B::Type{<:AbstractBond},::Type{Nothing})\noidtype(::Val{:Fock},B::Type{<:AbstractBond},::Type{<:Table})\n\nGet the compatible Fock OID type with an AbstractBond type and a Table/Nothing type.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.DegreesOfFreedom.otype-Tuple{Val{:Fock},Type{#s235} where #s235<:(Term{\'F\',Species,V,N,C,A,M} where M<:Union{Nothing, TermModulate} where A<:TermAmplitude where C<:TermCouplings where N where V<:Number where Species),Type{#s234} where #s234<:AbstractBond,Type{#s233} where #s233<:Union{Nothing, Dict{I,Int64} where I<:Index}}",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.DegreesOfFreedom.otype",
    "category": "method",
    "text": "otype(::Val{:Fock},O::Type{<:Term{\'F\'}},B::Type{<:AbstractBond},T::Type{<:Union{Nothing,Table}})\notype(::Val{:Fock},O::Type{<:Term{\'B\'}},B::Type{<:AbstractBond},T::Type{<:Union{Nothing,Table}})\n\nGet the compatible Fock operator type with a Term type, an AbstractBond type and a Table/Nothing type.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.isnormalordered-Tuple{FockOperator}",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.isnormalordered",
    "category": "method",
    "text": "isnormalordered(opt::FockOperator) -> Bool\n\nJudge whether a FockOperator is normal ordered.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.nambufockindextotuple",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.nambufockindextotuple",
    "category": "function",
    "text": "nambufockindextotuple\n\nIndicate that the filtered attributes are (:scope,:nambu,:site,:orbital,:spin) when converting a Fock index to tuple.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.usualfockindextotuple",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.usualfockindextotuple",
    "category": "function",
    "text": "usualfockindextotuple\n\nIndicate that the filtered attributes are (:scope,:site,:orbital,:spin) when converting a Fock index to tuple.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.σʸ-Tuple{String}",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.σʸ",
    "category": "method",
    "text": "σʸ(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing) -> Couplings{ID{<:NTuple{2,FCID}},FockCoupling{2,Int,ID{<:NTuple{2,FCID}}}}\n\nThe Pauli matrix σʸ, which can act on the space of spins(\"sp\"), orbitals(\"ob\"), sublattices(\"sl\") or particle-holes(\"ph\").\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.σˣ-Tuple{String}",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.σˣ",
    "category": "method",
    "text": "σˣ(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing) -> Couplings{ID{<:NTuple{2,FCID}},FockCoupling{2,Int,ID{<:NTuple{2,FCID}}}}\n\nThe Pauli matrix σˣ, which can act on the space of spins(\"sp\"), orbitals(\"ob\"), sublattices(\"sl\") or particle-holes(\"ph\").\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.σᶻ-Tuple{String}",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.σᶻ",
    "category": "method",
    "text": "σᶻ(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing) -> Couplings{ID{<:NTuple{2,FCID}},FockCoupling{2,Int,ID{<:NTuple{2,FCID}}}}\n\nThe Pauli matrix σᶻ, which can act on the space of spins(\"sp\"), orbitals(\"ob\"), sublattices(\"sl\") or particle-holes(\"ph\").\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.σ⁰-Tuple{String}",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.σ⁰",
    "category": "method",
    "text": "σ⁰(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing) -> Couplings{ID{<:NTuple{2,FCID}},FockCoupling{2,Int,ID{<:NTuple{2,FCID}}}}\n\nThe Pauli matrix σ⁰, which can act on the space of spins(\"sp\"), orbitals(\"ob\"), sublattices(\"sl\") or particle-holes(\"ph\").\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.σ⁺-Tuple{String}",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.σ⁺",
    "category": "method",
    "text": "σ⁺(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing) -> Couplings{ID{<:NTuple{2,FCID}},FockCoupling{2,Int,ID{<:NTuple{2,FCID}}}}\n\nThe Pauli matrix σ⁺, which can act on the space of spins(\"sp\"), orbitals(\"ob\"), sublattices(\"sl\") or particle-holes(\"ph\").\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.FockPackage.σ⁻-Tuple{String}",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.FockPackage.σ⁻",
    "category": "method",
    "text": "σ⁻(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing) -> Couplings{ID{<:NTuple{2,FCID}},FockCoupling{2,Int,ID{<:NTuple{2,FCID}}}}\n\nThe Pauli matrix σ⁻, which can act on the space of spins(\"sp\"), orbitals(\"ob\"), sublattices(\"sl\") or particle-holes(\"ph\").\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.Terms.statistics-Tuple{BOperator}",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.Terms.statistics",
    "category": "method",
    "text": "statistics(opt::BOperator)\nstatistics(::Type{<:BOperator})\n\nGet the statistics of BOperator.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Essentials.Terms.statistics-Tuple{FOperator}",
    "page": "Fock package",
    "title": "Hamiltonian.Essentials.Terms.statistics",
    "category": "method",
    "text": "statistics(opt::FOperator) -> Char\nstatistics(::Type{<:FOperator}) -> Char\n\nGet the statistics of FOperator.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Prerequisites.Interfaces.:⊗-Union{Tuple{N}, Tuple{FockCoupling{N,V,I,OS,SS} where SS<:Subscripts where OS<:Subscripts where I<:(ID{#s254} where #s254<:Tuple{Vararg{FCID,N}}) where V<:Number,FockCoupling{N,V,I,OS,SS} where SS<:Subscripts where OS<:Subscripts where I<:(ID{#s254} where #s254<:Tuple{Vararg{FCID,N}}) where V<:Number}} where N",
    "page": "Fock package",
    "title": "Hamiltonian.Prerequisites.Interfaces.:⊗",
    "category": "method",
    "text": "⊗(fc1::FockCoupling{N},fc2::FockCoupling{N}) where N -> FockCoupling\n\nGet the direct product between two Fock couplings.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Hamiltonian.Prerequisites.Interfaces.expand-Union{Tuple{S}, Tuple{FockCoupling,PID,Fock}, Tuple{FockCoupling,PID,Fock,Union{Nothing, Val{S}}}} where S",
    "page": "Fock package",
    "title": "Hamiltonian.Prerequisites.Interfaces.expand",
    "category": "method",
    "text": "expand(fc::FockCoupling,pid::PID,fock::Fock,species::Union{Val{S},Nothing}=nothing) where S -> Union{FCExpand,Tuple{}}\nexpand(fc::FockCoupling,pids::NTuple{R,PID},focks::NTuple{R,Fock},species::Union{Val{S},Nothing}=nothing) where {R,S} -> Union{FCExpand,Tuple{}}\n\nExpand a Fock coupling with the given set of point ids and Fock degrees of freedom.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Base.:*-Tuple{FockCoupling,FockCoupling}",
    "page": "Fock package",
    "title": "Base.:*",
    "category": "method",
    "text": "*(fc1::FockCoupling,fc2::FockCoupling) -> FockCoupling\n\nGet the multiplication between two Fock couplings.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Base.adjoint-Tuple{FID}",
    "page": "Fock package",
    "title": "Base.adjoint",
    "category": "method",
    "text": "adjoint(fid::FID) -> FID\n\nGet the adjoint of a Fock id.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Base.repr-Tuple{FockCoupling}",
    "page": "Fock package",
    "title": "Base.repr",
    "category": "method",
    "text": "repr(fc::FockCoupling) -> String\n\nGet the repr representation of a Fock coupling.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Base.show-Tuple{IO,FockCoupling}",
    "page": "Fock package",
    "title": "Base.show",
    "category": "method",
    "text": "show(io::IO,fc::FockCoupling)\n\nShow a Fock coupling.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Base.union-Union{Tuple{P}, Tuple{Type{P},Type{FID}}} where P<:PID",
    "page": "Fock package",
    "title": "Base.union",
    "category": "method",
    "text": "union(::Type{P},::Type{FID}) where P<:PID\n\nGet the union type of PID and FID.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/FockPackage.html#Manual-1",
    "page": "Fock package",
    "title": "Manual",
    "category": "section",
    "text": "Modules=[FockPackage]\nOrder=  [:module,:constant,:type,:macro,:function]"
},

{
    "location": "man/Essentials/SpinPackage.html#",
    "page": "Spin package",
    "title": "Spin package",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Essentials.SpinPackage"
},

{
    "location": "man/Essentials/SpinPackage.html#Spin-package-1",
    "page": "Spin package",
    "title": "Spin package",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/SpinPackage.html#Spin-degrees-of-freedom-1",
    "page": "Spin package",
    "title": "Spin degrees of freedom",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/SpinPackage.html#SID-and-Spin-1",
    "page": "Spin package",
    "title": "SID and Spin",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/SpinPackage.html#SIndex-1",
    "page": "Spin package",
    "title": "SIndex",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/SpinPackage.html#Predefined-spin-operators-1",
    "page": "Spin package",
    "title": "Predefined spin operators",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/SpinPackage.html#Spin-terms-1",
    "page": "Spin package",
    "title": "Spin terms",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/SpinPackage.html#SCID-and-SpinCoupling-1",
    "page": "Spin package",
    "title": "SCID and SpinCoupling",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/SpinPackage.html#Predefined-spin-couplings-1",
    "page": "Spin package",
    "title": "Predefined spin couplings",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/SpinPackage.html#Predefined-spin-terms-1",
    "page": "Spin package",
    "title": "Predefined spin terms",
    "category": "section",
    "text": ""
},

{
    "location": "man/Essentials/SpinPackage.html#Hamiltonian.Essentials.SpinPackage.SCID",
    "page": "Spin package",
    "title": "Hamiltonian.Essentials.SpinPackage.SCID",
    "category": "type",
    "text": "SCID(;center=wildcard,atom=wildcard,orbital=wildcard,tag=\'i\',subscript=wildcard)=SCID(center,atom,orbital,tag,subscript)\n\nThe id of a spin coupling.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/SpinPackage.html#Hamiltonian.Essentials.SpinPackage.SID",
    "page": "Spin package",
    "title": "Hamiltonian.Essentials.SpinPackage.SID",
    "category": "type",
    "text": "SID <: IID\n\nThe spin id.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/SpinPackage.html#Hamiltonian.Essentials.SpinPackage.SID-Tuple{}",
    "page": "Spin package",
    "title": "Hamiltonian.Essentials.SpinPackage.SID",
    "category": "method",
    "text": "SID(;orbital::Int=1,spin::Real=0.5,tag::Char=\'i\')\n\nCreate a spin id.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/SpinPackage.html#Hamiltonian.Essentials.SpinPackage.SIndex",
    "page": "Spin package",
    "title": "Hamiltonian.Essentials.SpinPackage.SIndex",
    "category": "type",
    "text": "SIndex{S} <: Index{PID{S},SID}\n\nThe spin index.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/SpinPackage.html#Hamiltonian.Essentials.SpinPackage.SOperator",
    "page": "Spin package",
    "title": "Hamiltonian.Essentials.SpinPackage.SOperator",
    "category": "type",
    "text": "SOperator(value::Number,id::ID{<:NTuple{N,OID}}) where N\n\nSpin operator.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/SpinPackage.html#Hamiltonian.Essentials.SpinPackage.Spin",
    "page": "Spin package",
    "title": "Hamiltonian.Essentials.SpinPackage.Spin",
    "category": "type",
    "text": "Spin <: Internal{SID}\n\nThe spin interanl degrees of freedom.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/SpinPackage.html#Hamiltonian.Essentials.SpinPackage.Spin-Tuple{}",
    "page": "Spin package",
    "title": "Hamiltonian.Essentials.SpinPackage.Spin",
    "category": "method",
    "text": "Spin(;atom::Int=1,norbital::Int=1,spin::Real=0.5)\n\nConstruct a spin degrees of freedom.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/SpinPackage.html#Hamiltonian.Essentials.SpinPackage.SpinCoupling",
    "page": "Spin package",
    "title": "Hamiltonian.Essentials.SpinPackage.SpinCoupling",
    "category": "type",
    "text": "SpinCoupling(value::Number,id::ID{<:NTuple{N,SCID}},subscripts::Subscripts) where N\nSpinCoupling{N}(    value::Number=1;\n                    tags::NTuple{N,Char},\n                    centers::Union{NTuple{N,Int},Nothing}=nothing,\n                    atoms::Union{NTuple{N,Int},Nothing}=nothing,\n                    orbitals::Union{NTuple{N,Int},Subscript,Nothing}=nothing\n                    ) where N\n\nSpin coupling.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/SpinPackage.html#Hamiltonian.Essentials.SpinPackage.SpinTerm",
    "page": "Spin package",
    "title": "Hamiltonian.Essentials.SpinPackage.SpinTerm",
    "category": "type",
    "text": "SpinTerm(   id::Symbol,value::Number;\n            neighbor::Int,\n            couplings::Union{Tuple{<:Tuple{Vararg{Couplings}},<:Function},Coupling,Couplings},\n            amplitude::Union{Function,Nothing}=nothing,\n            modulate::Union{Function,Bool}=false,\n            factor::Number=1\n            )\n\nSpin term.\n\nType alias for Term{\'B\',:SpinTerm,<:Number,Int,<:TermCouplings,<:TermAmplitude,<:Union{TermModulate,Nothing}}.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/SpinPackage.html#Hamiltonian.Essentials.DegreesOfFreedom.oidtype-Tuple{Val{:Spin},Type{#s237} where #s237<:AbstractBond,Type{Nothing}}",
    "page": "Spin package",
    "title": "Hamiltonian.Essentials.DegreesOfFreedom.oidtype",
    "category": "method",
    "text": "oidtype(::Val{:Spin},B::Type{<:AbstractBond},::Type{Nothing})\noidtype(::Val{:Spin},B::Type{<:AbstractBond},::Type{<:Table})\n\nGet the compatible spin OID type with an AbstractBond type and a Table/Nothing type.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/SpinPackage.html#Hamiltonian.Essentials.DegreesOfFreedom.otype-Tuple{Val{:Spin},Type{#s235} where #s235<:(Term{\'B\',Species,V,N,C,A,M} where M<:Union{Nothing, TermModulate} where A<:TermAmplitude where C<:TermCouplings where N where V<:Number where Species),Type{#s234} where #s234<:AbstractBond,Type{#s233} where #s233<:Union{Nothing, Dict{I,Int64} where I<:Index}}",
    "page": "Spin package",
    "title": "Hamiltonian.Essentials.DegreesOfFreedom.otype",
    "category": "method",
    "text": "otype(::Val{:Spin},O::Type{<:Term{\'B\'}},B::Type{<:AbstractBond},T::Type{<:Union{Nothing,Table}})\n\nGet the compatible spin operator type with a Term type, an AbstractBond type and a Table/Nothing type.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/SpinPackage.html#Hamiltonian.Essentials.SpinPackage.Heisenberg-Tuple{}",
    "page": "Spin package",
    "title": "Hamiltonian.Essentials.SpinPackage.Heisenberg",
    "category": "method",
    "text": "Heisenberg(;centers::Union{NTuple{2,Int},Nothing}=nothing,\n            atoms::Union{NTuple{2,Int},Nothing}=nothing,\n            orbitals::Union{NTuple{2,Int},Subscript,Nothing}=nothing\n            ) -> Couplings{ID{<:NTuple{2,SCID}},SpinCoupling{2,Float,ID{<:NTuple{2,SCID}}}}\n\nThe Heisenberg couplings.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/SpinPackage.html#Hamiltonian.Essentials.SpinPackage.Ising-Tuple{Char}",
    "page": "Spin package",
    "title": "Hamiltonian.Essentials.SpinPackage.Ising",
    "category": "method",
    "text": "Ising(  tag::Char;\n        centers::Union{NTuple{2,Int},Nothing}=nothing,\n        atoms::Union{NTuple{2,Int},Nothing}=nothing,\n        orbitals::Union{NTuple{2,Int},Subscript,Nothing}=nothing\n        ) -> Couplings{ID{<:NTuple{2,SCID}},SpinCoupling{2,Float,ID{<:NTuple{2,SCID}}}}\n\nThe Ising couplings.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/SpinPackage.html#Hamiltonian.Essentials.SpinPackage.Sʸ-Tuple{}",
    "page": "Spin package",
    "title": "Hamiltonian.Essentials.SpinPackage.Sʸ",
    "category": "method",
    "text": "Sʸ(;atom::Union{Int,Nothing}=nothing,orbital::Union{Int,Nothing}=nothing) -> Couplings{ID{<:NTuple{2,SCID}},SpinCoupling{2,Float,ID{<:NTuple{2,SCID}}}}\n\nThe single Sʸ coupling.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/SpinPackage.html#Hamiltonian.Essentials.SpinPackage.Sˣ-Tuple{}",
    "page": "Spin package",
    "title": "Hamiltonian.Essentials.SpinPackage.Sˣ",
    "category": "method",
    "text": "Sˣ(;atom::Union{Int,Nothing}=nothing,orbital::Union{Int,Nothing}=nothing) -> Couplings{ID{<:NTuple{2,SCID}},SpinCoupling{2,Float,ID{<:NTuple{2,SCID}}}}\n\nThe single Sˣ coupling.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/SpinPackage.html#Hamiltonian.Essentials.SpinPackage.Sᶻ-Tuple{}",
    "page": "Spin package",
    "title": "Hamiltonian.Essentials.SpinPackage.Sᶻ",
    "category": "method",
    "text": "Sᶻ(;atom::Union{Int,Nothing}=nothing,orbital::Union{Int,Nothing}=nothing) -> Couplings{ID{<:NTuple{2,SCID}},SpinCoupling{2,Float,ID{<:NTuple{2,SCID}}}}\n\nThe single Sᶻ coupling.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/SpinPackage.html#Hamiltonian.Essentials.SpinPackage.usualspinindextotuple",
    "page": "Spin package",
    "title": "Hamiltonian.Essentials.SpinPackage.usualspinindextotuple",
    "category": "function",
    "text": "usualspinindextotuple\n\nIndicate that the filtered attributes are (:scope,:site,:orbital) when converting a spin index to tuple.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/SpinPackage.html#Hamiltonian.Essentials.Terms.statistics-Tuple{SOperator}",
    "page": "Spin package",
    "title": "Hamiltonian.Essentials.Terms.statistics",
    "category": "method",
    "text": "statistics(opt::SOperator) -> Char\nstatistics(::Type{<:SOperator}) -> Char\n\nGet the statistics of SOperator.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/SpinPackage.html#Hamiltonian.Prerequisites.Interfaces.expand-Union{Tuple{S}, Tuple{SpinCoupling,PID,Spin}, Tuple{SpinCoupling,PID,Spin,Union{Nothing, Val{S}}}} where S",
    "page": "Spin package",
    "title": "Hamiltonian.Prerequisites.Interfaces.expand",
    "category": "method",
    "text": "expand(sc::SpinCoupling,pid::PID,spin::Spin,species::Union{Val{S},Nothing}=nothing) where S -> Union{SCExpand,Tuple{}}\nexpand(sc::SpinCoupling,pids::NTuple{N,PID},spins::NTuple{N,Spin},species::Union{Val{S},Nothing}=nothing) where {N,S} -> Union{SCExpand,Tuple{}}\n\nExpand a spin coupling with the given set of point ids and spin degrees of freedom.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/SpinPackage.html#Hamiltonian.Prerequisites.Interfaces.matrix",
    "page": "Spin package",
    "title": "Hamiltonian.Prerequisites.Interfaces.matrix",
    "category": "function",
    "text": "matrix(sid::SID,dtype::Type{<:Number}=Complex{Float}) -> Matrix{dtype}\n\nGet the matrix representation of a sid.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/SpinPackage.html#Base.:*-Tuple{SpinCoupling,SpinCoupling}",
    "page": "Spin package",
    "title": "Base.:*",
    "category": "method",
    "text": "*(sc1::SpinCoupling,sc2::SpinCoupling) -> SpinCoupling\n\nGet the multiplication between two spin couplings.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/SpinPackage.html#Base.adjoint-Tuple{SID}",
    "page": "Spin package",
    "title": "Base.adjoint",
    "category": "method",
    "text": "adjoint(sid::SID) -> SID\n\nGet the adjoint of a spin id.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/SpinPackage.html#Base.repr-Tuple{SpinCoupling}",
    "page": "Spin package",
    "title": "Base.repr",
    "category": "method",
    "text": "repr(sc::SpinCoupling) -> String\n\nGet the repr representation of a spin coupling.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/SpinPackage.html#Base.show-Tuple{IO,SpinCoupling}",
    "page": "Spin package",
    "title": "Base.show",
    "category": "method",
    "text": "show(io::IO,sc::SpinCoupling)\n\nShow a spin coupling.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/SpinPackage.html#Base.union-Union{Tuple{P}, Tuple{Type{P},Type{SID}}} where P<:PID",
    "page": "Spin package",
    "title": "Base.union",
    "category": "method",
    "text": "union(::Type{P},::Type{SID}) where {P<:PID}\n\nGet the union type of PID and SID.\n\n\n\n\n\n"
},

{
    "location": "man/Essentials/SpinPackage.html#Manual-1",
    "page": "Spin package",
    "title": "Manual",
    "category": "section",
    "text": "Modules=[SpinPackage]\nOrder=  [:module,:constant,:type,:macro,:function]"
},

]}
