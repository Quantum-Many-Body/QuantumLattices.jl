```@meta
CurrentModule=Hamiltonian
```

# Hamiltonian

*Julia package for constructing and solving the Hamiltonians of quantum lattice systems.*

We provide a general framework to construct the **symbolic representation** of the Hamiltonian of any **quantum lattice system**, with the inputs as simple as its description by natural language. Based on this symbolic representation, we implement several algorithms, such as **TBA**, **ED**, **CPT/VCA**, **DMRG**, etc., to solve the quantum lattice system. Generic interfaces are offered to access to these algorithms. Only minor modifications need be made when the user alters from an algorithm to another.

## Package Features

* **Unitcell Description Framework**: by telling the information of the quantum lattice system within a unitcell, the construction of the symbolic representation of the Hamiltonian is just as simple as describing the system in a usual research paper;
* **Generic Engine-App Interfaces**: by regarding the relation between algorithms and tasks as that between engines and apps, automatic project management is realized, including that of result recording, data caching, parameter updating, code logging, task dependency, etc, furthermore, all algorithms are initialized in quite similiar ways with only minor modifications needed.

## Supported Algorithms
* **TBA**: tight-binding approximation for fermionic/bosonic systems;
* **ED**: exact diagonalizaiton for fermionic/hard-core-bosonic/spin systems;
* **CPT/VCA**: cluster perturbation theory and variational cluster approach for fermionic systems;
* **DMRG**: density matrix renormalization group for fermionic/hard-core-bosonic/spin systems;
* **FBFM**: spin wave theory for flatband ferromagnets.

## Python counterpart
[HamiltonianPy](https://github.com/waltergu/HamiltonianPy): in fact, the authors of this Julia package worked on the python package at first and only turned to Julia later.
