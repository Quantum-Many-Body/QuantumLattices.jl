push!(LOAD_PATH,"../src/")

using Documenter
using QuantumLattices

makedocs(
    format=             :html,
    clean=              false,
    sitename=           "QuantumLattices.jl",
    pages=              [
                        "Home"      =>  "index.md",
                        "Tutorial"  =>  [
                                "tutorial/UnitcellDescription.md",
                                "tutorial/EngineAppInterface.md",
                                ],
                        "Manual"    =>  [
                                "man/Interfaces.md",
                                "Prerequisites" => [
                                        "man/Prerequisites/Introduction.md",
                                        "man/Prerequisites/TypeTraits.md",
                                        "man/Prerequisites/Factories.md",
                                        "man/Prerequisites/CompositeStructures.md",
                                        "man/Prerequisites/SimpleTrees.md",
                                        "man/Prerequisites/NamedVectors.md",
                                        ],
                                "Mathematics" => [
                                        "man/Mathematics/Introduction.md",
                                        "man/Mathematics/Combinatorics.md",
                                        "man/Mathematics/VectorSpaces.md",
                                        "man/Mathematics/AlgebraOverFields.md",
                                        "man/Mathematics/QuantumNumbers.md",
                                        ],
                                "Essentials" => [
                                        "man/Essentials/Introduction.md",
                                        "man/Essentials/Spatials.md",
                                        "man/Essentials/DegreesOfFreedom.md",
                                        "man/Essentials/Terms.md",
                                        "man/Essentials/FockPackage.md",
                                        "man/Essentials/SpinPackage.md",
                                        "man/Essentials/Extensions.md",
                                        ],
                                ],
                        ],
    html_canonical=     "https://quantum-many-body.github.io/QuantumLattices.jl/latest/",
    assets=             ["assets/favicon.ico"]
)

deploydocs(
    repo=       "github.com/Quantum-Many-Body/QuantumLattices.jl.git",
    target=     "build",
    julia=      "1.0",
    osname=     "linux",
    deps=       nothing,
    make=       nothing,
)
