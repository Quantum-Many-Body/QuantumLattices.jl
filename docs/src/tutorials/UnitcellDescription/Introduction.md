```@meta
CurrentModule = QuantumLattices
```

# Introduction

A quantum lattice system can be completely described within its unitcell. Basically, this description should contain three types of information:

1) the spatial information, such as the coordinates of the points contained in the unitcell;
2) the internal degrees of freedom, such as the local Hilbert space on each point;
3) the couplings among different degrees of freedom, such as the interaction terms in the Hamiltonian.

In theory, as long as the above information is told, one could easily write down the operator representation of the Hamiltonian of the system. For example, in the phrase *"the single orbital electronic Hubbard model with only nearest neighbor hopping on a one dimensional lattice with only two sites"*, *"one dimensional lattice with only two sites"* is the spatial information, *"single orbital electronic"* defines the local Hilbert spaces, and *"Hubbard model with only nearest neighbor hopping"* describes the terms present in the Hamiltonian. From this phrase, we also know that the Hamiltonian of the system is

```math
H=tc^†_{1↑}c_{2↑}+tc^†_{2↑}c_{1↑}+tc^†_{1↓}c_{2↓}+tc^†_{2↓}c_{1↓}+Uc^†_{1↑}c_{1↑}c^†_{1↓}c_{1↓}+Uc^†_{2↑}c_{2↑}c^†_{2↓}c_{2↓}
```

where ``t`` is the hopping amplitude， ``U`` is the Hubbard interaction strength and the electronic annihilation/creation operator carries a site index and a spin index. Actually, the **unitcell description framework** follows exactly after the above train of thought. For example, the aforementioned system can be constructed by the following codes

```@example
using QuantumLattices
using SymPy: symbols

# define the unitcell
lattice = Lattice("L2P", [Point(PID(1), (0.0,)), Point(PID(2), (1.0,))])

# define the internal degrees of freedom
config = IDFConfig{Fock}(pid->Fock(norbital=1, nspin=2, nnambu=2), lattice.pids)

# define the terms
t = Hopping{'F'}(:t, symbols("t", real=true), 1)
U = Hubbard{'F'}(:U, symbols("U", real=true))

# get the Hamiltonian
operators = expand(Generator((t, U), Bonds(lattice), config, nothing, false))
```
The last line displays all the generated operators in the Hamiltonian in the latex form. Here, in the subscript of the electronic annihilation/creation operator, an extra orbital index is also displayed. In the following sections listed below, we will explain in brief how these codes work. Specifically, in [Spatial info of a unitcell](@ref), we will introduce how to construct the unitcell as . For more detailed explanations, the manual of [`Essentials`](@ref essentials) can also be referred.

```@contents
Pages = [
        "SpatialInfoOfAUnitcell.md",
        "InternalDegreesOfFreedom.md",
        "CouplingsAmongDifferentDegreesOfFreedom.md",
        "GeneratorOfOperators.md",
        ]
Depth = 2
```
