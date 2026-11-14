# Quantum Gravity Context for the SU(2) Lattice Laboratory

## Purpose

This laboratory studies SU(2) (and a prototype SU(3)) lattice gauge theory with
the Wilson action. The same gauge group, the same link variables and the same
plaquette operators appear in a second field: loop quantum gravity (LQG), where
SU(2) is the structure group of the Ashtekar–Barbero connection.

This note records the dictionary between the two settings, so that results
computed in this repository can be read in that language, and says clearly
where the dictionary stops.

**Nothing in this repository is a model of quantum gravity.** The lab runs a
Euclidean Yang–Mills Monte Carlo on a fixed hypercubic lattice. This note does
not change that. It explains a correspondence of kinematics.

## The dictionary

| Lattice gauge theory (this lab) | Loop quantum gravity |
|---|---|
| Link variable \(U_\ell \in SU(2)\) | Holonomy \(h_e \in SU(2)\) of the connection along an edge |
| Plaquette \(U_p\) (product around a square) | Holonomy around a small loop, a discrete curvature \(F_{ab}\) |
| Wilson loop \(\operatorname{Tr}\,U_C\) | Loop state, the basic gauge-invariant observable |
| Gauge-invariant functions of links | Spin-network functions |
| Irreducible representation \(j\) on a link | Spin label \(j\) on an edge of a spin network |
| Electric field Casimir \(E^2 = j(j+1)\) (Kogut–Susskind) | Squared area flux through the surface the edge punctures |
| Gauss law at a vertex | Invariance of the intertwiner at a node |
| Wilson action \(\beta\sum_p(1-\tfrac12\operatorname{Re}\operatorname{Tr}U_p)\) | (Not part of the kinematics; a gauge-theory Hamiltonian) |

The Kogut–Susskind Hamiltonian for SU(2) is

\[
H
=
\frac{g^2}{2}\sum_{\ell} j_\ell(j_\ell+1)
+
\frac{1}{g^2}\sum_p\bigl(2-\operatorname{Tr}U_p\bigr),
\]

with \(\beta = 4/g^2\) in the convention used by this lab's Wilson action
\(\beta\sum_p(1-\tfrac12\operatorname{Re}\operatorname{Tr}U_p)\). The first term
is diagonal in the spin basis. The second multiplies by group elements and
changes spins. Both are already finite-dimensional once the spins are
truncated.

## The one number worth knowing

In LQG the area operator on a surface punctured by edges with spins \(j_i\)
has eigenvalues

\[
A = 8\pi\gamma\,\ell_P^2\sum_i\sqrt{j_i(j_i+1)},
\]

where \(\gamma\) is the Barbero–Immirzi parameter. The lattice gauge theory has
the same Casimir \(j(j+1)\) in its electric term. The two are the same
mathematical object read in different physical units. No coupling between
them is claimed.

For a single spin-\(\tfrac12\) puncture, \(\sqrt{j(j+1)}=\sqrt3/2\) and

\[
A_{1/2} = 4\sqrt3\,\pi\gamma\,\ell_P^2 .
\]

With the value \(\gamma=\ln 2/(\pi\sqrt3)\approx 0.1274\) obtained from
black-hole entropy counting with spin-\(\tfrac12\) punctures only, this is
\(4\ln 2\;\ell_P^2\approx 2.7726\,\ell_P^2\). The all-spin counting gives a
different value (about 0.2375, quoted from the literature, not rederived here).
The value of \(\gamma\) is model-dependent and is not predicted by this lab.

## What the lab can contribute

Because the lab already implements the pieces that appear on the lattice side,
a few numerical statements are available at no extra cost:

- **Link Casimir energies.** In the strong-coupling limit the electric term
  dominates, and low-lying states are labelled by spins on links. One link
  carrying spin \(\tfrac12\) costs \(\tfrac{g^2}{2}\cdot\tfrac34=\tfrac{3g^2}{8}\)
  in the convention above. A closed spin-\(\tfrac12\) loop around one
  plaquette (four links) therefore costs \(\tfrac{3g^2}{2}\), which is the
  strong-coupling gap quoted in the companion note `SU2_LATTICE_GAP.md` of the
  `openai_math_yang_mills` fork. That is a Hamiltonian statement; this lab's
  Monte Carlo does not compute it.
- **Exact two-dimensional results.** In two dimensions the Wilson area law is
  exact; see [`spin_networks_and_wilson_loops.md`](spin_networks_and_wilson_loops.md).
  This is a clean validation target for the Wilson-loop and Creutz-ratio code.
- **Gauge invariance and Gauss law checks.** The lab's gauge-transformation
  tests (`gauge_invariance.md`) are the lattice version of invariance of a
  spin-network state under local \(SU(2)\) rotations.

## What the dictionary does not give

- **No diffeomorphism invariance.** The lab uses a fixed hypercubic lattice
  with fixed spacing. LQG states live on graphs embedded up to diffeomorphism.
- **No dynamics of geometry.** The lattice is a background. In LQG the
  Hamiltonian constraint is the open problem; the Wilson-action Hamiltonian is
  not that constraint.
- **Euclidean sampling, not canonical quantization.** Monte Carlo samples the
  Euclidean path integral. LQG kinematics is canonical.
- **Different SU(2).** In the lab SU(2) is a Yang–Mills gauge group with a
  dimensionless coupling \(g\). In LQG it is the rotation group of the spatial
  triad, and \(\gamma\) plays the role the coupling does here. Weak isospin of
  the Standard Model is a third, unrelated SU(2).
- **No continuum, no gravity limit.** The lab makes finite-lattice
  statements only. See [`quantum_gravity_limits.md`](quantum_gravity_limits.md).

## Related reading in this repository

- [`lattice_gauge_theory.md`](lattice_gauge_theory.md): lattice setup and Wilson action.
- [`gauge_invariance.md`](gauge_invariance.md): gauge transformations and invariants.
- [`wilson_loops.md`](wilson_loops.md) and [`creutz_ratios.md`](creutz_ratios.md): observables.
- [`spin_networks_and_wilson_loops.md`](spin_networks_and_wilson_loops.md): gauge-invariant states and exact counts.
- [`quantum_gravity_limits.md`](quantum_gravity_limits.md): what is and is not claimed.
