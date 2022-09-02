using Documenter
using DocumenterTools
using QuantumLattices

makedocs(
    format=     Documenter.HTML(
                    prettyurls= !("local" in ARGS),
                    canonical=  "https://quantum-many-body.github.io/QuantumLattices.jl/latest/",
                    assets=     ["assets/favicon.ico"],
                    analytics=  "UA-89508993-1",
                ),
    sitename=   "QuantumLattices.jl",
    pages=      [
                    "Home"      =>  "index.md",
                    "Tutorials" =>  [
                        "Unitcell Description" => [
                            "tutorials/UnitcellDescription/Introduction.md",
                            "tutorials/UnitcellDescription/SpatialInfoOfAUnitcell.md",
                            "tutorials/UnitcellDescription/InternalDegreesOfFreedom.md",
                            "tutorials/UnitcellDescription/CouplingsAmongDifferentDegreesOfFreedom.md",
                            "tutorials/UnitcellDescription/GeneratorOfOperators.md",
                        ],
                        "Advanced Usage" => [
                            "tutorials/AdvancedUsage/Introduction.md",
                            "tutorials/AdvancedUsage/LaTeXFormatOutputs.md",
                            "tutorials/AdvancedUsage/BoundaryConditions.md",
                            "tutorials/AdvancedUsage/Transformations.md",
                            "tutorials/AdvancedUsage/ProjectManagement.md",
                        ]
                    ],
                    "Manual"    =>  [
                        "man/Interfaces.md",
                        "Prerequisites" => [
                            "man/Prerequisites/Introduction.md",
                            "man/Prerequisites/Combinatorics.md",
                            "man/Prerequisites/Traits.md",
                            "man/Prerequisites/CompositeStructures.md",
                            "man/Prerequisites/SimpleTrees.md",
                            "man/Prerequisites/NamedVectors.md",
                            "man/Prerequisites/VectorSpaces.md",
                        ],
                        "Essentials" => [
                            "man/Essentials/Introduction.md",
                            "man/Essentials/QuantumOperators.md",
                            "man/Essentials/QuantumNumbers.md",
                            "man/Essentials/Spatials.md",
                            "man/Essentials/DegreesOfFreedom.md",
                            "man/Essentials/QuantumSystems.md",
                            "man/Essentials/Frameworks.md",
                        ],
                    ]
                ]
)

deploydocs(
    repo=       "github.com/Quantum-Many-Body/QuantumLattices.jl.git",
    target=     "build",
    deps=       nothing,
    make=       nothing,
)
