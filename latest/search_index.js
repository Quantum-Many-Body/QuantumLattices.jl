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
    "text": "Utilities of the Hamiltonian package.Pages=  [\n        \"NamedVector.md\",\n        \"GoodQuantumNumber.md\",\n        ]\nDepth=2"
},

{
    "location": "man/Utilities/NamedVector.html#",
    "page": "NamedVector",
    "title": "NamedVector",
    "category": "page",
    "text": "CurrentModule=Hamiltonian.Utilities.NamedVector"
},

{
    "location": "man/Utilities/NamedVector.html#Hamiltonian.Utilities.NamedVector.AbstractNamedVector",
    "page": "NamedVector",
    "title": "Hamiltonian.Utilities.NamedVector.AbstractNamedVector",
    "category": "type",
    "text": "AbstractNamedVector{T,A<:AbstractVector{T}}\n\nAbstract type for all concrete named vectors. To subtype it, please note:\n\nThe concrete types must have the field values::A, which is used to store the values by design. A recommended template for the subtype is  struct YourNamedVector{T,A} <: AbstractNamedVector{T,A}      values::A  end\nThe concrete types must implement their own Base.fieldnames, which defines the names of the type. A recommended template for this method is  Base.fieldnames(::Type{YourNamedVector},private=false)=private ? tuple(YourNames...,:values) : YourNames\nArithmetic operations, such as +, -, *, \\, %, ÷, etc. are NOT supported. Thus users have to overload these operations themselves if needed.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Hamiltonian.Utilities.NamedVector.@namedvector",
    "page": "NamedVector",
    "title": "Hamiltonian.Utilities.NamedVector.@namedvector",
    "category": "macro",
    "text": "@namedvector typename fieldnames dtype::Union{Expr,Symbol}=:nothing atype::Union{Expr,Symbol}=:nothing supertypename=:AbstractNamedVector\n\nConstruct a concrete named vector with the type name being typename, fieldnames specified by fieldnames, and optionally, the type parameters specified by dtype and atype, and the supertype specified by supertypename.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.:==-Tuple{Hamiltonian.Utilities.NamedVector.AbstractNamedVector,Hamiltonian.Utilities.NamedVector.AbstractNamedVector}",
    "page": "NamedVector",
    "title": "Base.:==",
    "category": "method",
    "text": "==(nv1::AbstractNamedVector,nv2::AbstractNamedVector)\n\nOverloaded == operator.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.eltype-Union{Tuple{Type{NV}}, Tuple{NV}, Tuple{A}, Tuple{T}, Tuple{Type{NV},Int64}} where NV<:Hamiltonian.Utilities.NamedVector.AbstractNamedVector{T,A} where A where T",
    "page": "NamedVector",
    "title": "Base.eltype",
    "category": "method",
    "text": "eltype(::Type{NV},choice::Int=1) where NV<:AbstractNamedVector{T,A} where {T,A}\n\nGet the type parameter of a concrete AbstractNamedVector.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.getindex-Tuple{Hamiltonian.Utilities.NamedVector.AbstractNamedVector,Int64}",
    "page": "NamedVector",
    "title": "Base.getindex",
    "category": "method",
    "text": "getindex(nv::AbstractNamedVector,index::Int)\n\nGet the value by the [] syntax.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.getproperty-Tuple{Hamiltonian.Utilities.NamedVector.AbstractNamedVector,Symbol}",
    "page": "NamedVector",
    "title": "Base.getproperty",
    "category": "method",
    "text": "getproperty(nv::AbstractNamedVector,key::Symbol)\n\nGet the value by the . syntax.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.hash-Tuple{Hamiltonian.Utilities.NamedVector.AbstractNamedVector,UInt64}",
    "page": "NamedVector",
    "title": "Base.hash",
    "category": "method",
    "text": "hash(nv::AbstractNamedVector,h::UInt)\n\nHash a concrete AbstractNamedVector.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.iterate",
    "page": "NamedVector",
    "title": "Base.iterate",
    "category": "function",
    "text": "iterate(nv::AbstractNamedVector)\niterate(nv::AbstractNamedVector,state)\n\nIterate over the values of a concrete AbstractNamedVector.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.keys-Tuple{Hamiltonian.Utilities.NamedVector.AbstractNamedVector}",
    "page": "NamedVector",
    "title": "Base.keys",
    "category": "method",
    "text": "keys(nv::AbstractNamedVector)\n\nIterate over the names.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.length-Tuple{Hamiltonian.Utilities.NamedVector.AbstractNamedVector}",
    "page": "NamedVector",
    "title": "Base.length",
    "category": "method",
    "text": "length(nv::AbstractNamedVector)\n\nGet the length of a concrete AbstractNamedVector.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.pairs-Tuple{Hamiltonian.Utilities.NamedVector.AbstractNamedVector}",
    "page": "NamedVector",
    "title": "Base.pairs",
    "category": "method",
    "text": "pairs(nv::AbstractNamedVector)\n\nIterate over the name=>value pairs.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.replace-Tuple{Hamiltonian.Utilities.NamedVector.AbstractNamedVector}",
    "page": "NamedVector",
    "title": "Base.replace",
    "category": "method",
    "text": "replace(nv::AbstractNamedVector;kwargs...)\n\nReturn a copy of a concrete AbstractNamedVector with some of the filed values replaced by the keyword arguments.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.setindex!-Tuple{Hamiltonian.Utilities.NamedVector.AbstractNamedVector,Any,Int64}",
    "page": "NamedVector",
    "title": "Base.setindex!",
    "category": "method",
    "text": "setindex!(nv::AbstractNamedVector,value,index::Int)\n\nSet the value by the [] syntax.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.setproperty!-Tuple{Hamiltonian.Utilities.NamedVector.AbstractNamedVector,Symbol,Any}",
    "page": "NamedVector",
    "title": "Base.setproperty!",
    "category": "method",
    "text": "setproperty!(nv::AbstractNamedVector,key::Symbol,value)\n\nSet the value by the . syntax. Note the value of the field :values as a whole can not be changed directly by design. Instead, the elements of :values can be changed.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.show-Tuple{IO,Hamiltonian.Utilities.NamedVector.AbstractNamedVector}",
    "page": "NamedVector",
    "title": "Base.show",
    "category": "method",
    "text": "show(io::IO,nv::AbstractNamedVector)\n\nShow a concrete AbstractNamedVector.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.values-Tuple{Hamiltonian.Utilities.NamedVector.AbstractNamedVector}",
    "page": "NamedVector",
    "title": "Base.values",
    "category": "method",
    "text": "values(nv::AbstractNamedVector)\n\nIterate over the values.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#Base.zero-Tuple{Hamiltonian.Utilities.NamedVector.AbstractNamedVector}",
    "page": "NamedVector",
    "title": "Base.zero",
    "category": "method",
    "text": "zero(nv::AbstractNamedVector)\n\nGet a concrete AbstractNamedVector with all values being zero.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/NamedVector.html#NamedVector-1",
    "page": "NamedVector",
    "title": "NamedVector",
    "category": "section",
    "text": "A named vector is similiar to a named tuple, which associate each of its values with a name, represented with a Symbol. Unlike named tuples, the values of a named vector can be modified. Yet the names of a named vector cannot be changed. With experiences of namedtuple in Python, we treat the names of a named vector as type traits, by overloading the Base.filednames function, so that all instances of a certain concrete named vector share the same names, which is different from the predefined NamedTuple in Julia. For users who want to define their own concrete named vector, we recommend the use of the macro @namedvector.Modules=[NamedVector]\nOrder=  [:module,:constant,:type,:macro,:function]"
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
    "text": "The abstract type for the complete set of independent good quantum numbers for a single basis.Main features include:function fieldnames: get the names of the quantum numbers\nfunction periods: get the periods of the quantum numbers\narithmetic operations: +,-,*,⊕\nhashable: concrete instances can be used as keys for a dict or a set\niterable: concrete instances are iterable over their valuesIn particular, QuantumNumber <: AbstractNamedVector{Float64}, all features supported by AbstractNamedVector are also available for QuantumNumber. See also AbstractNamedVector.For convenience, 4 kinds of good quantum numbers are predefined in this module, i.e.SQN: for spin z-component reserved systems\nPQN: for particle number reserved systems\nSPQN: for both particle number and spin-z component reserved systems\nZ2QN: for systems with a Z_2 conservation quantum numberUsers who want to define their own Z_N-like quantum numbers must handle the periodicities in the construction function, otherwise, wrong results will be get when arithmetic operations, such as + or -, are involved. It is highly recommended to use the macro @quantumnumber to define your own concrete QuantumNumbers."
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#QuantumNumbers-1",
    "page": "Good quantum numbers",
    "title": "QuantumNumbers",
    "category": "section",
    "text": "The whole quantum numbers for the total bases.To achieve high efficiency:The quantum numbers are stored in a compressed form similiar to that of a CSC/CSR sparse matrix.\nThe contents of a QuantumNumbers are not an array of some kind of concrete QuantumNumbers, but an 2d array of floats with the columns being the values of those concrete QuantumNumbers.Main features include:function eltype: get the concrete type of the quantum numbers it contains\narithmetic operations: +,-,*,^,⊗,⊕\niteration: concrete QuantumNumbers it contains will be constructed and returned"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.PQN",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.PQN",
    "category": "type",
    "text": "PQN(N::Real)\n\nThe concrete QuantumNumber of a quantum system with particle number N conserved.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.QNSProtocol",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.QNSProtocol",
    "category": "type",
    "text": "Protocol used for the initilization of QuantumNumbers.\n\nqnscounts: initilization by counts\nqnsindptr: initilization by indptr\n\n\n\n\n\n"
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
    "text": "QuantumNumbers(form::Char,::Type{QN},contents::Array{Float64,2},info::Array{Int,1},protocol::QNSProtocol=qnscounts,regularize::Bool=true) where QN<:QuantumNumber\nQuantumNumbers(form::Char,contents::Array{QN,1},info::Array{Int,1},protocol::QNSProtocol=qnscounts) where QN<:QuantumNumber\nQuantumNumbers(qn::QuantumNumber,count::Int=1)\n\nThe whole quantum numbers of the total bases of a Hilbert space.\n\n\n\n\n\n"
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
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.:⊕",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.:⊕",
    "category": "function",
    "text": "⊕(qns::QN...;signs::Union{Array{Int,1},Nothing}=nothing) where QN<:QuantumNumber\n⊕(qnses::QuantumNumbers{QN}...;signs::Union{Array{Int,1},Nothing}=nothing) where QN\n\nGet the direct sum of some QuantumNumbers or QuantumNumberss.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.:⊗",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.:⊗",
    "category": "function",
    "text": "⊗(::Type{QN},qn1::QuantumNumber,qn2::QuantumNumber) where QN<:QuantumNumber\n⊗(qnses::QuantumNumbers{QN}...;signs::Union{Array{Int,1},Nothing}=nothing) where QN\n\nGet the direct product of some QuantumNumbers or QuantumNumberss.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.regularization-Union{Tuple{N}, Tuple{QN}, Tuple{Type{QN},AbstractArray{Float64,N}}} where N where QN<:QuantumNumber",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.regularization",
    "category": "method",
    "text": "regularization(::Type{QN},array::AbstractArray{Float64,N}) where {QN<:QuantumNumber,N}\n\nGet the regularized array of the array representation of a QuantumNumber/QuantumNumbers by the periods.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Hamiltonian.Utilities.GoodQuantumNumber.regularize!",
    "page": "Good quantum numbers",
    "title": "Hamiltonian.Utilities.GoodQuantumNumber.regularize!",
    "category": "function",
    "text": "regularize!(::Type{QN},array::AbstractVector{Float64}) where QN<:QuantumNumber\nregularize!(::Type{QN},array::AbstractArray{Float64,2}) where QN<:QuantumNumber\n\nRegularize an array representation of a QuantumNumber/QuantumNumbers by the concrete QuantumNumber\'s periods.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.:*",
    "page": "Good quantum numbers",
    "title": "Base.:*",
    "category": "function",
    "text": "*(qn::QuantumNumber,factor::Integer)\n*(factor::Integer,qn::QuantumNumber)\n*(qns::QuantumNumbers,factor::Integer)\n*(factor::Integer,qns::QuantumNumbers)\n\nOverloaded * operator for QuantumNumber and QuantumNumbers..\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.:+",
    "page": "Good quantum numbers",
    "title": "Base.:+",
    "category": "function",
    "text": "+(qn::QuantumNumber)\n+(qn1::QN,qn2::QN) where QN<:QuantumNumber\n+(qn1::QN,qn2::QN,qns::QN...) where QN<:QuantumNumber\n+(qns::QuantumNumbers)\n+(qn::QN,qns::QuantumNumbers{QN}) where QN\n+(qns::QuantumNumbers{QN},qn::QN) where QN\n\nOverloaded + operator for QuantumNumber and QuantumNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.:-",
    "page": "Good quantum numbers",
    "title": "Base.:-",
    "category": "function",
    "text": "-(qn::QuantumNumber)\n-(qn1::QN,qn2::QN) where QN<:QuantumNumber\n-(qns::QuantumNumbers)\n-(qn::QN,qns::QuantumNumbers{QN}) where QN\n-(qns::QuantumNumbers{QN},qn::QN) where QN\n\nOverloaded - operator for QuantumNumber and QuantumNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.:^",
    "page": "Good quantum numbers",
    "title": "Base.:^",
    "category": "function",
    "text": "^(qn::QuantumNumber,factor::Integer)\n^(qns::QuantumNumbers,factor::Integer)\n\nOverloaded ^ operator for QuantumNumber and QuantumNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.eltype-Union{Tuple{Type{#s32} where #s32<:QuantumNumbers{QN}}, Tuple{QN}} where QN",
    "page": "Good quantum numbers",
    "title": "Base.eltype",
    "category": "method",
    "text": "eltype(::Type{<:QuantumNumbers{QN}}) where QN\n\nGet the type of the concrete QuantumNumber contained in QuantumNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.getindex-Tuple{QuantumNumbers,Int64}",
    "page": "Good quantum numbers",
    "title": "Base.getindex",
    "category": "method",
    "text": "getindex(qns::QuantumNumbers,index::Int)\n\nOverloaded [] operator.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.iterate",
    "page": "Good quantum numbers",
    "title": "Base.iterate",
    "category": "function",
    "text": "iterate(qns::QuantumNumbers,state::Int=1)\n\nIterate over the concrete QuantumNumbers the QuantumNumbers contains.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.keys-Tuple{QuantumNumbers}",
    "page": "Good quantum numbers",
    "title": "Base.keys",
    "category": "method",
    "text": "keys(qns::QuantumNumbers)\n\nIterate over the concrete QuantumNumbers the QuantumNumbers contains.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.length",
    "page": "Good quantum numbers",
    "title": "Base.length",
    "category": "function",
    "text": "length(qns::QuantumNumbers,duplicate::Bool=true)\n\nGet the number of qunatum numbers in the QuantumNumbers.\n\nduplicate==true: the duplicate quantum numbers are counted duplicately. Then the result equals the dimension of the QuantumNumbers.\nduplicate==false: only unduplicate quantum numbers are counted. Then the result equals the number of columns of the QuantumNumbers\' contents.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.pairs",
    "page": "Good quantum numbers",
    "title": "Base.pairs",
    "category": "function",
    "text": "pairs(qns::QuantumNumbers,protocol::QNSProtocol=qnsindptr)\n\nIterate over the QuantumNumber=>slice or QuantumNumber=>count pairs.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.show-Tuple{IO,QuantumNumbers}",
    "page": "Good quantum numbers",
    "title": "Base.show",
    "category": "method",
    "text": "show(io::IO,qns::QuantumNumbers)\n\nShow a QuantumNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.string-Tuple{QuantumNumbers}",
    "page": "Good quantum numbers",
    "title": "Base.string",
    "category": "method",
    "text": "string(qns::QuantumNumbers)\n\nConvert a QuantumNumbers to string.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Base.values",
    "page": "Good quantum numbers",
    "title": "Base.values",
    "category": "function",
    "text": "values(qns::QuantumNumbers,protocol::QNSProtocol=qnsindptr)\n\nIterate over the slices/counts of the QuantumNumbers.\n\n\n\n\n\n"
},

{
    "location": "man/Utilities/GoodQuantumNumber.html#Manual-1",
    "page": "Good quantum numbers",
    "title": "Manual",
    "category": "section",
    "text": "Modules=[GoodQuantumNumber]\nOrder=  [:module,:constant,:type,:macro,:function]"
},

]}
