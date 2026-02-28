```@meta
CurrentModule = QuantumLattices
```

# [Introduction](@id UnitcellDescriptionIntroduction)

A quantum lattice system can be completely described based on its unitcell. Essentially, this description requires three types of information:

1. The **spatial information**, such as the coordinates of the points contained in the unitcell;
2. The **internal degrees of freedom**, such as the local algebra acting on the local Hilbert space at each point;
3. The **couplings among different degrees of freedom**, such as the terms present in the Hamiltonian.

In theory, given the above information, one can easily write down the operator representation of the Hamiltonian of the system. For example, in the phrase *"the single orbital electronic Hubbard model with only nearest neighbor hopping on a one dimensional lattice with only two sites"*, *"one dimensional lattice with only two sites"* provides the spatial information, *"single orbital electronic"* defines the local Hilbert space and thus the local algebra, and *"Hubbard model with only nearest neighbor hopping"* describes the terms present in the Hamiltonian. From this description, we can derive that the Hamiltonian of the system is

```math
H=tc^\dagger_{1\uparrow}c_{2\uparrow}+tc^\dagger_{2\uparrow}c_{1\uparrow}+tc^\dagger_{1\downarrow}c_{2\downarrow}+tc^\dagger_{2\downarrow}c_{1\downarrow}+Uc^\dagger_{1\uparrow}c_{1\uparrow}c^\dagger_{1\downarrow}c_{1\downarrow}+Uc^\dagger_{2\uparrow}c_{2\uparrow}c^\dagger_{2\downarrow}c_{2\downarrow}
```

where ``t`` is the hopping amplitude, ``U`` is the Hubbard interaction strength, and the electronic creation/annihilation operator $c^\dagger_{i\sigma}/c_{i\sigma}$ carries a site index $i$ ($i=1, 2$) and a spin index $\sigma$ ($\sigma=\uparrow, \downarrow$). The **unitcell description framework** follows exactly this train of thought. For example, the aforementioned system can be constructed with the following code:

```@example
using QuantumLattices
using SymPy: Sym, symbols

# define the unitcell
lattice = Lattice([zero(Sym)], [one(Sym)])

# define the internal degrees of freedom, i.e., the single-orbital spin-1/2 fermionic algebra
hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(lattice))

# define the terms
t = Hopping(:t, symbols("t", real=true), 1)
U = Hubbard(:U, symbols("U", real=true))

# get the Hamiltonian
operators = expand(OperatorGenerator(bonds(lattice, 1), hilbert, (t, U)))
```
The last line displays all the generated operators in the Hamiltonian in LaTeX format. Here, in the subscript of the electronic annihilation/creation operator, an extra orbital index is also displayed.

In the following pages, we will explain in detail how these codes work:
* [Spatial information of a unitcell](@ref): Introduces the construction of the unitcell and the ways to obtain bonds at different orders of nearest neighbors.
* [Internal degrees of freedom](@ref): Explains the hierarchy of the internal degrees of freedom and how they are organized for different categories of quantum systems.
* [Couplings among different degrees of freedom](@ref): Discusses the ways to specify the terms present in the Hamiltonian.
* [Generator of operators](@ref): Shows how to combine everything to get the operator representation of the Hamiltonian of a quantum lattice system.

For more detailed information, please refer to the manual of this package.

```@contents
Pages = [
        "SpatialInfoOfAUnitcell.md",
        "InternalDegreesOfFreedom.md",
        "CouplingsAmongDifferentDegreesOfFreedom.md",
        "GeneratorOfOperators.md",
        ]
Depth = 2
```
