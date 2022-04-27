using QuantumLattices

lattice = Lattice(:L2P, [Point(PID(1), [0.0]), Point(PID(2), [1.0])])
hilbert = Hilbert(pid=>Fock{:f}(norbital=1, nspin=2) for pid in lattice.pids)
t = Hopping(:t, 1.0, 1)
U = Hubbard(:U, 8.0)
operators = expand(Generator((t, U), Bonds(lattice), hilbert))