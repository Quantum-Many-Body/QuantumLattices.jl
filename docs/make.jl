push!(LOAD_PATH,"../src/")

using Documenter
using Hamiltonian

makedocs(
    format=             :html,
    clean=              false,
    sitename=           "Hamiltonian.jl",
    pages=              [
                        "Home"      =>  "index.md",
                        "Tutorial"  =>  [
                                "tutorial/UnitcellDescription.md",
                                "tutorial/EngineAppInterface.md",
                                ],
                        "Manual"    =>  [
                                "Utilities" =>  [
                                        "man/Utilities/Introduction.md",
                                        "man/Utilities/Interface.md",
                                        "man/Utilities/TypeTrait.md",
                                        "man/Utilities/Factory.md",
                                        "man/Utilities/CompositeStructure.md",
                                        "man/Utilities/Tree.md",
                                        "man/Utilities/NamedVector.md",
                                        "man/Utilities/Combinatorics.md",
                                        "man/Utilities/AlgebraOverField.md",
                                        "man/Utilities/QuantumNumber.md",
                                        ],
                                "Essentials" => [
                                        "man/Essentials/Introduction.md",
                                        "man/Essentials/Spatial.md",
                                        "man/Essentials/DegreeOfFreedom.md",
                                        "man/Essentials/FockPackage.md",
                                        ],
                                ],
                        ],
    html_canonical=     "https://quantum-many-body.github.io/Hamiltonian.jl/latest/",
    assets=             ["assets/favicon.ico"]
)

deploydocs(
    repo=       "github.com/Quantum-Many-Body/Hamiltonian.jl.git",
    target=     "build",
    julia=      "1.0",
    osname=     "linux",
    deps=       nothing,
    make=       nothing,
)
