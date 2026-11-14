# Spin Networks, Wilson Loops and Exact Counts

## Purpose

[`quantum_gravity_context.md`](quantum_gravity_context.md) gives the dictionary
between lattice gauge theory and loop quantum gravity kinematics. This note
works out three concrete, checkable statements that sit on that dictionary and
connect to observables this lab already computes:

1. how many gauge-invariant states a vertex supports,
2. what the two-dimensional Wilson area law is exactly,
3. why the Wilson "area law" and the LQG "area spectrum" must not be confused.

Every number below was recomputed numerically while writing this note (Haar
integration with `scipy`, exact diagonalisation with `numpy`). The scripts were
scratch checks and are not part of the repository.

## 1. Counting gauge-invariant states at a vertex

Gauss's law says the states at a vertex are invariant under the local gauge
group. For SU(2), the number of invariant vectors in a tensor product of
representations is a Haar integral of the product of characters:

\[
\dim\operatorname{Inv}\Bigl(\bigotimes_i V_{j_i}\Bigr)
=
\int_{SU(2)} dU\;\prod_i \chi_{j_i}(U).
\]

For SU(2), with \(U\) conjugate to \(\operatorname{diag}(e^{i\theta},e^{-i\theta})\),
the Haar measure on class functions is \(\tfrac{2}{\pi}\sin^2\theta\,d\theta\)
on \([0,\pi]\) and \(\chi_{1/2}=2\cos\theta\).

For \(2n\) spin-\(\tfrac12\) legs the count is the Catalan number
\(C_n=\tfrac{1}{n+1}\binom{2n}{n}\):

| legs \(2n\) | 2 | 4 | 6 | 8 | 10 | 12 | 14 |
|---|---|---|---|---|---|---|---|
| invariants | 1 | 2 | 5 | 14 | 42 | 132 | 429 |

Numerical integration reproduces all seven values to about \(10^{-6}\).

The first row with more than one state is the four-valent vertex. Its
two-dimensional intertwiner space is the smallest space on which a
non-trivial "volume" can be defined. For the operator
\(Q=\epsilon_{abc}J^{(1)}_aJ^{(2)}_bJ^{(3)}_c\) with \(J=\sigma/2\) on each leg,
exact diagonalisation on the invariant subspace of four spin-\(\tfrac12\) legs
gives eigenvalues \(\pm\sqrt3/4\approx\pm0.4330\). This is the "quantum
tetrahedron". The spectrum of \(Q\) is symmetric about zero, so the sign is an
orientation, and only \(|Q|\) is a volume-like quantity.

**Relevance to the lab.** The lab's links carry group elements, not spins, so
it never constructs these states. The counts tell you how large the
gauge-invariant Hilbert space is on a given graph if you work in the spin basis,
which is the setting of the exact Hamiltonian computations in the
`openai_math_yang_mills` fork (`SU2_LATTICE_GAP.md`). They do not enter the
Monte Carlo.

## 2. The exact area law in two dimensions

In two dimensions there is no transverse direction, plaquettes are independent
after gauge fixing, and the plaquette weight \(e^{(\beta/2)\operatorname{Tr}U}\)
(equivalent to the lab's \(e^{-S_W}\), up to a constant) expands in characters.
The ratio of the spin-\(\tfrac12\) coefficient to the trivial one is

\[
u(\beta)=\frac{I_2(\beta)}{I_1(\beta)},
\]

where \(I_\nu\) are modified Bessel functions. A rectangular Wilson loop of
area \(A\) in the fundamental representation then behaves as \(u(\beta)^{A}\) on
the infinite lattice, so the string tension is exactly

\[
\sigma(\beta)=-\ln\frac{I_2(\beta)}{I_1(\beta)}
=\frac{3}{2\beta}+\frac{3}{4\beta^2}+\frac{3}{16\beta^3}-\frac{9}{16\beta^4}+O(\beta^{-5}).
\]

Numerical check of the closed form against direct Haar integration of the
spin-\(\tfrac12\) character coefficient:

| \(\beta\) | direct integration | \(-\ln(I_2/I_1)\) | 4-term expansion |
|---|---|---|---|
| 0.5 | 2.08975 | 2.08975 | (not valid, \(\beta\) small) |
| 2.0 | 0.83672 | 0.83672 | 0.9258 |
| 5.0 | 0.32942 | 0.32942 | 0.3306 |

The first two columns agree to machine precision at every \(\beta\) tested. The
asymptotic series is only useful at large \(\beta\).

**Use as a validation target.** The Creutz ratio
\(\chi(R,R)=-\ln\frac{W(R,R)W(R-1,R-1)}{W(R,R-1)^2}\)
removes perimeter terms, and in two dimensions on an infinite lattice it equals
\(\sigma(\beta)\) for every \(R\). The lab's `creutz.py` and the dimension-general
SU(2) code (`dimensional_scaling.md`) can therefore be run at \(d=2\) and
compared with the closed form above, up to finite-volume and statistical error.
**This comparison has not been run.** It is the natural first experiment if
you want a quantitative test of the Wilson-loop pipeline, with a known answer
and no continuum assumption.

On a periodic two-dimensional lattice with \(A\) plaquettes the partition
function is, in character form,

\[
Z=\sum_{j}\Bigl(\frac{2\,I_{2j+1}(\beta)}{\beta}\Bigr)^{A},
\]

as derived in `SU2_LATTICE_GAP.md` of the `openai_math_yang_mills` fork (not
rederived here). The \(j=0\) term dominates for large \(A\), and the
\(j=\tfrac12\) term is suppressed by \(u(\beta)^{A}\), so plaquette-type
quantities have corrections that are exponentially small in the lattice area.
Finite-size effects for individual Wilson loops on a torus are not worked out
here.

## 3. Two different "areas"

| | Wilson area law | LQG area spectrum |
|---|---|---|
| Quantity | \(\ln W\propto -\sigma\cdot\text{Area}\) for large loops | \(A=8\pi\gamma\ell_P^2\sum\sqrt{j(j+1)}\) |
| Meaning of "area" | Minimal surface bounded by the loop, in lattice units | Surface in the 3-geometry punctured by spin-network edges |
| Coefficient | String tension \(\sigma\) (energy per length), dimensionful | \(\gamma\ell_P^2\), Planck units |
| Probes | Confinement of colour flux | Quantum geometry |

The two share only the Casimir-like structure \(j(j+1)\) in the strong-coupling
regime, where the area law coefficient comes from the lowest non-trivial
character. Using one as evidence for the other would be a category error.
The same warning applies to the entanglement-entropy "area law" that appears in
some of the preprints cited in
[`openai_preprint_references.md`](openai_preprint_references.md).

## 4. What would make this note stronger

- Run the \(d=2\) Creutz-ratio comparison above and record agreement or the
  size of the discrepancy.
- Extend the character-coefficient check to spin-1 (adjoint) loops. In two
  dimensions the same argument gives \(I_3(\beta)/I_1(\beta)\) per plaquette.
  In \(d\ge3\) the large-loop behaviour of adjoint loops is different, because
  gluons can screen the adjoint charge (a standard result, not tested here).
- State the same checks for SU(3), where the Weyl-integral counting of
  invariants needs the SU(3) Haar measure rather than the one-angle formula
  used here.

## Limitations

These are kinematic and exactly solvable statements. They do not address
four-dimensional dynamics, a continuum limit, or the existence of a mass gap,
and they do not bear on the Yang–Mills Millennium Problem. The 2D theory has
no gap in the infinite-volume sense.
