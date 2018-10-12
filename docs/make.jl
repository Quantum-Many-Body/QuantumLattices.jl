push!(LOAD_PATH,"../src/")

using Documenter
using Hamiltonian,Hamiltonian.Utilities.GoodQuantumNumber

makedocs(
    format=             :html,
    clean=              false,
    sitename=           "Hamiltonian.jl",
    modules=            [Hamiltonian,Hamiltonian.Utilities,Hamiltonian.Utilities.GoodQuantumNumber],
    pages=              [
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
                        ],
    html_canonical=     "https://quantum-many-body.github.io/Hamiltonian.jl/latest/"
)

# deploydocs(
#     deps=       Deps.pip("pygments","mkdocs","python-markdown-math"),
#     repo=       "github.com/Quantum-Many-Body/Hamiltonian.jl.git",
#     julia=      "1.0",
#     osname=     "linux"
# )

deploydocs(
    repo=       "github.com/Quantum-Many-Body/Hamiltonian.jl.git",
    target=     "build",
    julia=      "1.0",
    osname=     "linux",
    deps=       nothing,
    make=       nothing,
)
