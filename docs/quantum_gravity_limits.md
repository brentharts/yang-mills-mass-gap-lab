# Quantum Gravity: Scope, Limits and Possible Next Steps

## Purpose

This note fixes what this repository does and does not claim in connection with
quantum gravity, so that the two companion notes
([`quantum_gravity_context.md`](quantum_gravity_context.md) and
[`spin_networks_and_wilson_loops.md`](spin_networks_and_wilson_loops.md)) are
not read as more than they are.

## What is claimed

- SU(2) link variables, plaquettes, Wilson loops and the Casimir \(j(j+1)\)
  appear both in this lab and in loop quantum gravity kinematics, with the
  dictionary in `quantum_gravity_context.md`.
- Some exact finite statements on that dictionary are numerically
  reproducible: Catalan-number counts of vertex invariants, the
  four-valent triple-product eigenvalues \(\pm\sqrt3/4\), and the
  two-dimensional string tension \(\sigma(\beta)=-\ln(I_2/I_1)\).
- The lab's Wilson-loop and Creutz-ratio code could be validated against the
  two-dimensional closed form. That comparison is proposed, not done.

## What is not claimed

| Topic | Status |
|---|---|
| This lab is a model of quantum gravity | **No.** Fixed background lattice, Yang–Mills dynamics. |
| Diffeomorphism invariance | Absent. |
| A Hamiltonian constraint | Not constructed. The Wilson-action Hamiltonian is not it. |
| Volume operator beyond one four-valent vertex | Not built here. |
| Problem of time | Not addressed. |
| Matter coupled to spin networks | None. |
| Value of the Barbero–Immirzi parameter | Not computed in this lab. Three counting conventions are derived in note G13 of `brentharts/quantum-gravity` (0.1274, 0.2375, 0.2741); the value is a property of the counting. |
| Any link to the Yang–Mills mass gap problem | None. The notes concern kinematics and the exactly solvable \(d=2\) theory, which has no mass gap in the infinite-volume sense. |
| Continuum limit | Not taken. |

## Wording to use and to avoid

Safe:

- "The Casimir \(j(j+1)\) that appears in the Kogut–Susskind electric term is
  the same function that gives the LQG area eigenvalues."
- "In two dimensions the Wilson area law is exact, with
  \(\sigma=-\ln(I_2/I_1)\)."
- "The number of SU(2)-invariant states at a vertex of \(2n\) spin-\(\tfrac12\)
  legs is the Catalan number \(C_n\)."

Avoid:

- Any sentence in which the Wilson area law, the string tension or a lattice
  mass gap is presented as supporting or testing loop quantum gravity.
- Treating the lattice coupling \(g\) as the Immirzi parameter, or the lattice
  spacing as the Planck length.
- Calling a lattice Hamiltonian "the Hamiltonian constraint".
- Reading any black-hole-entropy value of \(\gamma\) (0.1274, 0.2375 or
  0.2741, depending on the count) as a prediction of this lab.

## Related exact models outside this repository

The `openai_math_yang_mills` fork contains a separate set of exact
finite-dimensional models (spin-basis Hamiltonians for the two-plaquette ladder,
a three-plaquette strip, an L-shaped cluster with a four-valent vertex, and
checks of Penrose-type inequalities) with their own documentation in
`QUANTUM_GRAVITY.md`. Those are computed there, not here, and I have not re-run
them for this note. If one of them is brought into this lab, it should come
with tests in `tests/` and an entry in `docs/validation.md` like the existing
observables.

## Possible next steps, in order of cost

1. **\(d=2\) string-tension check.** Run the existing dimension-general SU(2)
   Monte Carlo with \(d=2\), compute Creutz ratios, and compare with
   \(-\ln(I_2/I_1)\) at several \(\beta\). Cheap, uses existing code, and gives
   a quantitative validation of the Wilson-loop pipeline.
2. **Spin-basis exact diagonalisation on a small graph.** Reproduce the
   single-plaquette strong-coupling gap \(3g^2/2\) in this repository's own
   conventions. A new module and tests would be needed.
3. **Vertex-invariant counting for SU(3).** Needs the SU(3) Haar measure on
   the maximal torus; relevant to the existing SU(3) prototype.
4. **Anything involving a volume operator or a constraint.** Substantial new
   work with no connection to the current Monte Carlo; not recommended inside
   this lab.

## Literature

Cited from memory or from the fork's notes; check each entry before using it in
a paper.

- K. G. Wilson, "Confinement of quarks", *Phys. Rev. D* 10 (1974) 2445.
- J. Kogut and L. Susskind, "Hamiltonian formulation of Wilson's lattice gauge
  theories", *Phys. Rev. D* 11 (1975) 395.
- C. Rovelli and L. Smolin, "Discreteness of area and volume in quantum
  gravity", *Nucl. Phys. B* 442 (1995) 593.
- T. Thiemann, "Quantum spin dynamics (QSD)", gr-qc/9606089.
- M. Domagała and J. Lewandowski, "Black hole entropy from quantum geometry",
  and K. Meissner, "Black hole entropy in loop quantum gravity" (2004), for the
  two-state all-spin value of \(\gamma\); A. Ghosh and P. Mitra, gr-qc/0411035,
  for the all-projection count (the arXiv number is as given in note G13 of
  `brentharts/quantum-gravity`; not checked here).

For rigorous background on two-dimensional asymptotically free models, see
[`openai_preprint_references.md`](openai_preprint_references.md).
