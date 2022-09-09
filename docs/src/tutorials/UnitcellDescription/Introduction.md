```@meta
CurrentModule = QuantumLattices
```

# [Introduction](@id UnitcellDescriptionIntroduction)

A quantum lattice system can be completely described based on its unitcell. Basically, this description should contain three types of information:

1) the spatial information, such as the coordinates of the points contained in the unitcell;
2) the internal degrees of freedom, such as the local algebra acting on the local Hilbert space at each point;
3) the couplings among different degrees of freedom, such as the interaction terms in the Hamiltonian.

In theory, as long as the above information is told, one could easily write down the operator representation of the Hamiltonian of the system. For example, in the phrase *"the single orbital electronic Hubbard model with only nearest neighbor hopping on a one dimensional lattice with only two sites"*, *"one dimensional lattice with only two sites"* is the spatial information, *"single orbital electronic"* defines the local Hilbert space and thus the local algebra, and *"Hubbard model with only nearest neighbor hopping"* describes the terms present in the Hamiltonian. From this phrase, we also know that the Hamiltonian of the system is

```math
H=tc^†_{1↑}c_{2↑}+tc^†_{2↑}c_{1↑}+tc^†_{1↓}c_{2↓}+tc^†_{2↓}c_{1↓}+Uc^†_{1↑}c_{1↑}c^†_{1↓}c_{1↓}+Uc^†_{2↑}c_{2↑}c^†_{2↓}c_{2↓}
```

where ``t`` is the hopping amplitude， ``U`` is the Hubbard interaction strength and the electronic annihilation/creation operator $c^\dagger_{i\sigma}/c_{i\sigma}$ carries a site index $i$ ($i=1, 2$) and a spin index $\sigma$ ($\sigma=\uparrow, \downarrow$). Actually, the **unitcell description framework** follows exactly after the above train of thought. For example, the aforementioned system can be constructed by the following codes:

```@example
using QuantumLattices
using SymPy: Sym, symbols

# define the unitcell
lattice = Lattice([zero(Sym)], [one(Sym)])

# define the internal degrees of freedom, single-orbital spin-1/2 fermionic algebra
hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(lattice))

# define the terms
t = Hopping(:t, symbols("t", real=true), 1)
U = Hubbard(:U, symbols("U", real=true))

# get the Hamiltonian
operators = expand(OperatorGenerator((t, U), bonds(lattice, 1), hilbert))
```
The last line displays all the generated operators in the Hamiltonian in the latex form. Here, in the subscript of the electronic annihilation/creation operator, an extra orbital index is also displayed.

In the following sections listed below, we will explain in brief how these codes work. Firstly, in the section of [Spatial information of a unitcell](@ref), we will introduce the construction of the unitcell as well as the ways to obtain the bonds at different orders of nearest neighbors. Secondly, in the section of [Internal degrees of freedom](@ref), we will explain the hierarchy of the internal degrees of freedom and discuss in detail how they are organized for different categories of quantum systems. Thirdly, in the section of [Couplings among different degrees of freedom](@ref), we will discuss the ways to specify the terms present in the Hamiltonian. Finally, in the section of [Generator of operators](@ref), we will show how to combine all to get the operator representation of the Hamiltonian of a quantum lattice system. For more detailed explanations, the manual of [`Essentials`](@ref essentials) can also be referred.

```@contents
Pages = [
        "SpatialInfoOfAUnitcell.md",
        "InternalDegreesOfFreedom.md",
        "CouplingsAmongDifferentDegreesOfFreedom.md",
        "GeneratorOfOperators.md",
        ]
Depth = 2
```
