module QuantumLattices

import LinearAlgebra: ishermitian, rank, mul!

# interfaces
export dimension, id, ishermitian, rank, value
export ⊕, ⊗, add!, div!, mul!, sub!
export decompose, decompose!, expand, expand!, permute
export dtype, kind, reset!, update, update!
function id end
function value end
function dimension end
function ⊕ end
function ⊗ end
function add! end
function sub! end
function div! end
function expand end
function expand! end
function decompose end
function decompose! end
function permute end
function dtype end
function kind end
function update end
function update! end
function reset! end

# Toolkit
include("Toolkit.jl")
using .Toolkit

# QuantumOperators
include("QuantumOperators.jl")
using .QuantumOperators
export ID, LaTeX, Operator, OperatorPack, OperatorProd, Operators, OperatorSet, OperatorSum, OperatorUnit, QuantumOperator, Representation
export LinearFunction, LinearTransformation, Matrixization, Permutation, RankFilter, TabledUnitSubstitution, Transformation, UnitSubstitution
export idtype, ishermitian, isscalartype, latexname, latexformat, matrix, script, sequence

# QuantumNumbers
include("QuantumNumbers.jl")
using .QuantumNumbers
export AbelianNumber, AbelianNumbers, Momenta, Momentum, Momentum₁, Momentum₂, Momentum₃, ParticleNumber, QuantumNumber, SpinfulParticle, Sz
export findindex, periods, @abeliannumber

# Spatials
include("Spatials.jl")
using .Spatials
export azimuth, azimuthd, distance, isintratriangle, isonline, isparallel, issubordinate, interlinks, minimumlengths, polar, polard, reciprocals, rotate, translate, tile, volume
export AbstractLattice, Bond, BrillouinZone, Lattice, Neighbors, Point, ReciprocalSpace, ReciprocalZone, ReciprocalPath, bonds, bonds!, icoordinate, isintracell, rcoordinate, save, selectpath, shrink, ReciprocalCurve
export azimuth, azimuthd, direction, distance, isintratriangle, isonline, isparallel, issubordinate, interlinks, minimumlengths, polar, polard, reciprocals, rotate, translate, tile, volume
export AbstractLattice, Bond, BrillouinZone, Lattice, Neighbors, Point, ReciprocalSpace, ReciprocalZone, ReciprocalPath, bonds, bonds!, icoordinate, isintracell, rcoordinate, save, selectpath, shrink
export @hexagon_str, @line_str, @rectangle_str

# DegreesOfFreedom
include("DegreesOfFreedom.jl")
using .DegreesOfFreedom
export Boundary, CompositeIID, CompositeIndex, CompositeInternal, Component, Constraint, Coupling, Hilbert, IID, IIDSpace, Index, Internal, MatrixCoupling, Metric, OperatorUnitToTuple, SimpleIID, SimpleInternal, Table, Term
export plain, statistics, @indexes

# QuantumSystems
include("QuantumSystems.jl")
using .QuantumSystems
export annihilation, creation, latexofbosons, latexoffermions, latexofparticles, @σ_str, @L_str
export Coulomb, FID, Fock, FockTerm, Hopping, Hubbard, InterOrbitalInterSpin, InterOrbitalIntraSpin, Onsite, PairHopping, Pairing, SpinFlip, isannihilation, iscreation, isnormalordered
export latexofspins, SID, Spin, totalspin, @Γ_str, @Γ′_str, @DM_str, @Heisenberg_str, @Ising_str
export DM, Heisenberg, Ising, Kitaev, SingleIonAnisotropy, SpinTerm, Zeeman, Γ, Γ′
export latexofphonons, Elastic, PID, Phonon, Kinetic, Hooke, PhononTerm

# Frameworks
include("Frameworks.jl")
using .Frameworks
export eager, lazy, Action, Algorithm, AnalyticalExpression, Assignment, CategorizedGenerator, CompositeGenerator, Eager, ExpansionStyle, Frontend, Generator, Image, Lazy, OperatorGenerator, Parameters, initialize, prepare!, run!, save

end
