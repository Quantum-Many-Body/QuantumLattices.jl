push!(LOAD_PATH,"../src/")

using Documenter
using Hamiltonian,Hamiltonian.Utilities.GoodQuantumNumber

makedocs(
    format=     :html,
    sitename=   "Hamiltonian.jl",
    modules=    [Hamiltonian,Hamiltonian.Utilities,Hamiltonian.Utilities.GoodQuantumNumber],
    clean=      false,
    pages=      [
                "Home"      =>  "index.md",
                "Tutorial"  =>  [
                        "tutorial/Unitcell Description.md",
                        "tutorial/Engine App Interface.md",
                        ],
                "Manual"    =>  [
                        "Utilities" =>  [
                                "man/Utilities/Introduction.md",
                                "man/Utilities/NamedVector.md",
                                "man/Utilities/GoodQuantumNumber.md",
                                ],
                        ],
                ]
)
