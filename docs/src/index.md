```@meta
CurrentModule=QuantumLattices
```

# QuantumLattices

*Julia package for the construction of quantum lattice systems.*

We provide a general framework to construct the **symbolic representation** of the Hamiltonian of any **quantum lattice system**, with the inputs as simple as its description by natural language. This symbolic representation serves as the prerequisite of solving quantum many-body problems, based on which several algorithms, such as **TBA**(tight-bind approximation), **SCMF**(self-consistent mean field theory), **ED**(exact diagonalizaiton), **CPT/VCA**(cluster perturbation theory and variational cluster approach ), **DMRG**(density matrix renormalization group), etc. can be implemented. Generic interfaces are defined to give a unified access to these algorithms although their real implementations come in seperate packages. Only minor modifications need be made when the user alters from one algorithm to another.

## Package Features

* **Unitcell Description Framework**: by telling the information of the quantum lattice system within a unitcell, the construction of the symbolic representation of the Hamiltonian is just as simple as describing the system in a usual research paper;
* **Generic Engine-App Interfaces**: by regarding the relation between algorithms and tasks as that between engines and apps, automatic project management is realized, including that of result recording, data caching, parameter updating, code logging, task dependency, etc, furthermore, all algorithms are initialized in quite similiar ways with only minor modifications needed.

## Supported Algorithms

Concrete algorithms are implemented in seperate packages (still in progess):
* **TBA**: tight-binding approximation for fermionic/bosonic systems;
* **SCMF**: self-consistent mean field theory for fermionic systems;
* **ED**: exact diagonalizaiton for fermionic/hard-core-bosonic/spin systems;
* **CPT/VCA**: cluster perturbation theory and variational cluster approach for fermionic systems;
* **DMRG**: density matrix renormalization group for fermionic/hard-core-bosonic/spin systems;
* **FBFMSW**: spin wave theory for flatband ferromagnets.

## Python counterpart
[HamiltonianPy](https://github.com/waltergu/HamiltonianPy): in fact, the authors of this Julia package worked on the python package at first and only turned to Julia later.
