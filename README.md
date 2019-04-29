# QuantumLattices.jl

[![Build Status](https://api.travis-ci.org/Quantum-Many-Body/QuantumLattices.jl.svg?branch=master)](https://travis-ci.org/Quantum-Many-Body/QuantumLattices.jl)
[![codecov](https://codecov.io/gh/Quantum-Many-Body/QuantumLattices.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Quantum-Many-Body/QuantumLattices.jl)
[![][docs-latest-img]][docs-latest-url]
[![][docs-stable-img]][docs-stable-url]
[![996.icu](https://img.shields.io/badge/link-996.icu-red.svg)](https://996.icu)
[![LICENSE](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![LICENSE](https://img.shields.io/badge/license-Anti%20996-blue.svg)](https://github.com/996icu/996.ICU/blob/master/LICENSE)

Julia package for the construction of quantum lattice systems.

## Installation

In Julia **v1.1+**, please type `]` in the REPL to use the package mode, then type this command:

```julia
pkg> add QuantumLattices
```

## Introduction

The core of the package is the constructuion of the **symbolic representations** of **lattice Hamiltonians**. This is based on the following  mathematical observations:
* the operators contained in the lattice Hamiltonian act on **local Hilbert spaces**, and
* the operators contained in the lattice Hamiltonian form **algebras over the complex field**.

The first observation is the starting point of our construction framework (i.e. [**the unitcell description framewrok**](https://quantum-many-body.github.io/QuantumLattices.jl/dev/tutorial/UnitcellDescription/)) and the second is the mathmatical foundation of our implementation of the **symbolic computation**.

It is noted that our implementation of the symbolic computation only involves
1) the mathematical operations between a scalar and an operator, and
2) the mathematical operations between two operators.

The symbolic operations between two scalars are **not** implemented becase:
* in condensed matter physics, for many cases, only the numerical values of operators are important because the analytical expressions can be too complicated to analyze or they may even not exist;
* our construction process of the operators and their mathematical operations are **completely compatible with the [`SymPy`](https://github.com/JuliaPy/SymPy.jl) package**, therefore, a fully symbolic computation can be acheived by a simple combination of both.

Specically, we provide several general functions to deal with three common kinds of systems in condensed matter physics
* **canonical fermionic systems**
* **canonical bosonic systems**
* **SU(2) spin systems**

Thees functions can also be applied to other systems if you extend our "protocols", which won't cost you much effort.

This package can serve as the **prerequisites of all quantum lattice algorithms/approaches**, e.g. TBA, SCMF, ED, CPT/VCA, DMRG, etc. Some of these algorithms have being under development by the same authors of this package elsewhere.

For tutorials and munuals of this pacakge, please refer to its docs.

## Documentation
- [**LATEST**][docs-latest-url] &mdash; **documentation of the latest version.**
- [**STABLE**][docs-stable-url] &mdash; **documentation of the stable version.**


## Contact
waltergu1989@gmail.com


[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://quantum-many-body.github.io/QuantumLattices.jl/latest/
[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://quantum-many-body.github.io/QuantumLattices.jl/stable/
