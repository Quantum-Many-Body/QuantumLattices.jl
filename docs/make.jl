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
                                        "man/Utilities/Interfaces.md",
                                        "man/Utilities/TypeTraits.md",
                                        "man/Utilities/Factories.md",
                                        "man/Utilities/CompositeStructures.md",
                                        "man/Utilities/Trees.md",
                                        "man/Utilities/NamedVectors.md",
                                        "man/Utilities/Combinatorics.md",
                                        "man/Utilities/AlgebraOverFields.md",
                                        "man/Utilities/QuantumNumbers.md",
                                        ],
                                "Essentials" => [
                                        "man/Essentials/Introduction.md",
                                        "man/Essentials/Spatials.md",
                                        "man/Essentials/DegreesOfFreedom.md",
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
