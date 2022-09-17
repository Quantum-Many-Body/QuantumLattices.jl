using Documenter
using DocumenterTools
using QuantumLattices

makedocs(
    format=     Documenter.HTML(
                    prettyurls= get(ENV, "CI", "false") == "true",
                    canonical=  "https://quantum-many-body.github.io/QuantumLattices.jl/latest/",
                    assets=     ["assets/favicon.ico"],
                    analytics=  "UA-89508993-1",
                ),
    sitename=   "QuantumLattices.jl",
    pages=      [
                    "Home"      =>  "index.md",
                    "Unitcell Description Framework" => [
                        "unitcell description framework/Introduction.md",
                        "unitcell description framework/SpatialInfoOfAUnitcell.md",
                        "unitcell description framework/InternalDegreesOfFreedom.md",
                        "unitcell description framework/CouplingsAmongDifferentDegreesOfFreedom.md",
                        "unitcell description framework/GeneratorOfOperators.md",
                    ],
                    "Advanced Topics" => [
                        "advanced topics/Introduction.md",
                        "advanced topics/LaTeXFormatOutputs.md",
                        "advanced topics/HybridSystems.md",
                        "advanced topics/BoundaryConditions.md",
                        "advanced topics/Transformations.md",
                        "advanced topics/OrderIndexes.md",
                        "advanced topics/ManageProjects.md",
                    ],
                    "Manual"    =>  [
                        "man/Toolkit.md",
                        "man/QuantumOperators.md",
                        "man/QuantumNumbers.md",
                        "man/Spatials.md",
                        "man/DegreesOfFreedom.md",
                        "man/QuantumSystems.md",
                        "man/Frameworks.md",
                    ]
                ]
)

deploydocs(
    repo=       "github.com/Quantum-Many-Body/QuantumLattices.jl.git",
    target=     "build",
    deps=       nothing,
    make=       nothing,
)
