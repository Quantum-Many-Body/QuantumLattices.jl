using Documenter
using QuantumLattices

makedocs(
    format=     Documenter.HTML(
                    prettyurls = get(ENV, "CI", "false") == "true",
                    canonical = "https://quantum-many-body.github.io/QuantumLattices.jl/latest/",
                    assets = ["assets/favicon.ico", "assets/custom.css"],
                    analytics = "UA-89508993-1",
                    size_threshold_warn = 204800,
                ),
    sitename=   "QuantumLattices.jl",
    pages=      [
                    "Home" => "index.md",
                    "Tutorials" => [
                        "tutorials/1-introduction.md",
                        "tutorials/2-lattice-and-spatial.md",
                        "tutorials/3-internal-degrees-of-freedom.md",
                        "tutorials/4-operators.md",
                        "tutorials/5-couplings-and-terms.md",
                        "tutorials/6-latticemodel.md",
                        "tutorials/7-algorithm-interface.md",
                    ],
                    "Advanced Topics" => [
                        "advanced topics/Introduction.md",
                        "advanced topics/HybridSystems.md",
                        "advanced topics/BoundaryConditions.md",
                        "advanced topics/LinearTransformations.md",
                    ],
                    "Manual" => [
                        "man/Toolkit.md",
                        "man/QuantumOperators.md",
                        "man/Spatials.md",
                        "man/DegreesOfFreedom.md",
                        "man/QuantumSystems.md",
                        "man/Frameworks.md",
                    ],
                ]
)

deploydocs(
    repo=       "github.com/Quantum-Many-Body/QuantumLattices.jl.git",
    target=     "build",
    deps=       nothing,
    make=       nothing,
)
