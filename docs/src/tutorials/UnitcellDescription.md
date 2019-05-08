```@meta
CurrentModule=QuantumLattices
```

```@setup unitcell
push!(LOAD_PATH,"../../../src/")
using QuantumLattices
```

# Unitcell Description

A quantum lattice system can be completely described within its unitcell. Bascially, this description should contain three types of information:

1) the spatial information, such as the the coordinates of the points contained in the unitcell;
2) the internal degrees of freedom, such as the local Hilbert space on each point;
3) the couplings among different degrees of freedom, such as the interaction terms in the Hamiltonian.

In theory, as long as the above information is told, one could easily write down the operator representation of the Hamiltonian of the system. For example, in the phrase *"the single orbital electronic Hubbard model with only nearest neighbor hopping on a one dimensioanl lattice with only two sites"*, *"one dimensioanl lattice with only two sites"* is the spatial information, *"single orbital electronic"* defines the local Hilbert spaces, and *"Hubbard model with only nearest neighbor hopping"* describes the terms present in the Hamiltonian. From this phrase, we also know that the Hamiltonian of the system is

```math
H=tc^†_{1↑}c_{2↑}+tc^†_{2↑}c_{1↑}+tc^†_{1↓}c_{2↓}+tc^†_{2↓}c_{1↓}+Uc^†_{1↑}c_{1↑}c^†_{1↓}c_{1↓}+Uc^†_{2↑}c_{2↑}c^†_{2↓}c_{2↓}
```

where ``t`` is the hopping amplitude and ``U`` is the Hubbard interaction strength. Actually, the **unitcell description framework** follows exactly after the above train of thought. For example, the forementioned system can be constructed by the following codes

```@example
using QuantumLattices
using SymPy: symbols

# define the unitcell
lattice=Lattice("L2P",[Point(PID(1),(0.0,)),Point(PID(2),(1.0,))])

# define the internal degrees of freedom
config=IDFConfig{Fock}(pid->Fock(norbital=1,nspin=2,nnambu=2),lattice.pids)

# define the terms
t,U=symbols("t,U",real=true)
t=Hopping{'F'}(:t,t,1)
U=Hubbard{'F'}(:U,U)

# get the Hamiltonian
operators=expand(Generator((t,U),Bonds(lattice),config,nothing,false))
```
The last line displays all the generated operators in the Hamiltonian in the latex form. In the following sections, we will explain in brief how these codes work. For detailed explanations, please refer to the manual of [`Essentials`](@ref essentails).

## Spatial info of a unitcell


## Internal degrees of freedom


## Couplings among different degrees of freedom
