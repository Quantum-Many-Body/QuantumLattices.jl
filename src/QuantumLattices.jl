module QuantumLattices

include("Interfaces.jl")
include("Prerequisites/Prerequisites.jl")
include("Essentials/Essentials.jl")

using .Interfaces
using .Prerequisites
using .Essentials

# Interfaces
export id, value, rank, dimension
export ⊕, ⊗, ⋅, add!, sub!, mul!, div!
export expand, expand!, decompose, decompose!, permute

# Essentials
export dtype, kind, update!, reset!

# Essentials.QuantumAlgebras
export ID, idtype, sequence

# Essentials.QuantumNumbers
export AbelianNumber, AbelianNumbers, @abeliannumber, periods

# Essentials.Spatials
export distance, azimuthd, azimuth, polard, polar, volume, isparallel, isonline, isintratriangle, issubordinate, reciprocals, translate, rotate
export PID, CPID, Point, Bond, Lattice, SuperLattice, Cylinder, Bonds, BrillouinZone
export allbonds, zerothbonds, insidebonds, acrossbonds, intrabonds, interbonds
export pidtype, rcoord, icoord, isintracell, bonds!, bonds

# Essentials.DegreesOfFreedom
export CompositeIID, CompositeInternal, Hilbert, Index, OID, Operator, Operators, OIDToTuple, Table, LaTeX, Boundary
export iidtype, isHermitian, latexformat, twist, plain

# Essentials.Terms
export IIDSpace, Subscript, IIDConstrain, Coupling, Couplings, Term, Parameters, Generator
export abbr, ismodulatable, otype, @subscript_str, @couplings

# Essentials.Frameworks
export App, Engine, Assignment, Algorithm
export prepare!, register!, run!, dependences, rundependences!

# Essentials.QuantumSystems
## Canonical fermionic/bosonic systems
export ANNIHILATION, CREATION, MAJORANA, fdefaultlatex, bdefaultlatex, usualfockindextotuple, nambufockindextotuple
export FID, Fock, FockCoupling, Onsite, Hopping, Pairing, Hubbard, InterOrbitalInterSpin, InterOrbitalIntraSpin, SpinFlip, PairHopping, Coulomb
export statistics, isnormalordered, @σ⁰_str, @σˣ_str, @σʸ_str, @σᶻ_str, @σ⁺_str, @σ⁻_str, @fc_str

## SU(2) spin systems
export sdefaultlatex, usualspinindextotuple
export SID, Spin, SpinCoupling, SpinTerm, totalspin
export @heisenberg_str, @ising_str, @gamma_str, @dm_str, @sˣ_str, @sʸ_str, @sᶻ_str, @sc_str

## Phononic systems
export pndefaultlatex, usualphononindextotuple
export PNID, Phonon, PhononCoupling, PhononKinetic, PhononPotential
export @kinetic_str, @potential_str

end
