module QuantumLattices

import LinearAlgebra: ishermitian, rank, mul!, ⋅

# interfaces
include("interfaces.jl")
export dimension, id, ishermitian, rank, value
export ⊕, ⊗, ⋅, add!, div!, mul!, sub!
export decompose, decompose!, expand, expand!, permute
export dtype, kind, reset!, update, update!

# Toolkit
include("Toolkit.jl")
using .Toolkit

# QuantumOperators
include("QuantumOperators.jl")
using .QuantumOperators
export ID, Operator, OperatorPack, OperatorProd, Operators, OperatorSum, OperatorUnit, LaTeX, QuantumOperator
export AbstractSubstitution, AbstractUnitSubstitution, Identity, MatrixRepresentation, Numericalization, Permutation, RankFilter, Transformation, UnitSubstitution
export ishermitian, latexname, latexformat, matrix, matrix!, script, sequence

# QuantumNumbers
include("QuantumNumbers.jl")
using .QuantumNumbers
export AbelianNumber, AbelianNumbers, Momentum, Momentum₁, Momentum₂, Momentum₃, ParticleNumber, SpinfulParticle, SpinZ
export particlenumbers, periods, spinfulparticles, spinzs, @abeliannumber

# Spatials
include("Spatials.jl")
using .Spatials
export azimuth, azimuthd, distance, isintratriangle, isonline, isparallel, issubordinate, interlinks, minimumlengths, polar, polard, reciprocals, rotate, translate, tile, volume
export AbstractLattice, Bond, BrillouinZone, Lattice, Neighbors, Point, ReciprocalSpace, ReciprocalZone, ReciprocalPath, Segment, Translations, bonds!, bonds, icoordinate, isintracell, rcoordinate, @translations_str
export @hexagon_str, @line_str, @rectangle_str

# DegreesOfFreedom
include("DegreesOfFreedom.jl")
using .DegreesOfFreedom
export Boundary, CompositeIID, CompositeIndex, CompositeInternal, Component, Constraint, Coupling, Hilbert, IID, IIDSpace, Index, Internal, MatrixCoupling, Metric, OperatorUnitToTuple, SimpleIID, SimpleInternal, Table, Term
export plain, statistics, @indexes

# QuantumSystems
include("QuantumSystems.jl")
using .QuantumSystems
export σ⁰, σˣ, σʸ, σᶻ, σ⁺, σ⁻, annihilation, creation, latexofbosons, latexoffermions, latexofparticles
export Coulomb, FID, Fock, FockTerm, Hopping, Hubbard, InterOrbitalInterSpin, InterOrbitalIntraSpin, Onsite, PairHopping, Pairing, SpinFlip, isnormalordered
export latexofspins, SID, Spin, SpinTerm, totalspin, @dm_str, @gamma_str, @heisenberg_str, @ising_str
export latexofphonons, Elastic, PID, Phonon, Kinetic, Hooke, PhononTerm

# Frameworks
include("Frameworks.jl")
using .Frameworks
export Action, Algorithm, AnalyticalExpression, Assignment, CompositeGenerator, Entry, Frontend, Image, OperatorGenerator, Parameters, RepresentationGenerator, prepare!, run!, rundependences!, save

end
