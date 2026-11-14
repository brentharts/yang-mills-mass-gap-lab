# References from the OpenAI Math Release Preprints

This page lists preprints from the OpenAI Math Release (as collected in the
[openai_math_yang_mills](https://github.com/brentharts/openai_math_yang_mills)
fork) that can be cited as **background literature** for this laboratory. The
BibTeX entries are in [`openai_preprints.bib`](openai_preprints.bib), copied
verbatim from each preprint's own README.

## Read this first

**None of these papers is about four-dimensional Yang–Mills, lattice gauge
theory, Wilson loops, glueballs or overrelaxation.** A keyword scan of all 751
preprint PDFs found no lattice-gauge-theory paper. Every entry below is an
*analogy* or a *methodological* reference, and each entry says what it does
not give you. In particular:

- Nothing here supports a claim that this lab says anything about the
  continuum Yang–Mills mass gap. The lab's scope stays finite-lattice.
- The preprints are release documents, not peer-reviewed publications. I have
  read abstracts and introductions, not verified the proofs.
- Cite them as "OpenAI Math Release preprint, 2026", as the BibTeX does.

## A. Exactly solved mass gaps in two-dimensional asymptotically free models

Relevant to: [`mass_gap_estimation.md`](mass_gap_estimation.md),
[`beta_finite_size_scaling.md`](beta_finite_size_scaling.md),
[`spectroscopy.md`](spectroscopy.md).

| Key | What it proves | Where to cite |
|---|---|---|
| `OAI:Exact-mass-asymptotics-for-the-two-dimensional-O4-lattice-model-October-5-2026` | Exact low-temperature asymptotics of the lattice mass of the 2D O(4) model, defined as the gap of the full Osterwalder–Schrader transfer operator (including rotation-invariant local-observable sectors), `m(β) ~ c √β e^{-πβ}`. | `beta_finite_size_scaling.md`: as an example that a lattice mass gap in an asymptotically free model closes exponentially in the inverse coupling, so a fixed-`β` plateau is not a continuum statement. |
| `OAI:Sharp-mass-bounds-for-the-two-dimensional-O4-model-September-23-2026` | Two-sided bounds of the same `√β e^{-πβ}` form for large `β`, and a positive full gap at every finite `β > 0`. | `mass_gap_estimation.md`, "Scientific Limitation": a positive finite-lattice gap is the expected situation and is not evidence of a continuum gap. |
| `OAI:The-canonical-massive-continuum-limit-of-the-two-dimensional-O3-model-October-4-2026` | A continuum limit of the 2D O(3) model satisfying the Osterwalder–Schrader axioms, with a positive Hamiltonian gap. | `mass_gap_estimation.md`: the lab's "Purpose" paragraph. This is what a *completed* version of the lab's finite-lattice programme looks like in two dimensions. |
| `OAI:An-Isolated-Particle-Pole-for-the-Two-Dimensional-O3-Spin-Field-October-4-2026` | The lowest mass is an isolated atom of the spectral measure, separated by a gap from the rest. | `spectroscopy.md`: motivation for why the lab separates the ground-state mass from the continuum threshold in effective-mass plateaus. |
| `OAI:Exponential-decay-in-two-dimensional-classical-On-models-September-23-2026` | Exponential decay of two-point functions of the O(n) model for every `n ≥ 3` and every finite temperature. | `mass_gap_estimation.md`, "Correlation Decay": a rigorous instance of the exponential-decay ansatz the lab fits, in a different model. |

**Caveats.** These are O(N) sigma models (spins on a sphere), not gauge
theories. Asymptotic freedom is shared; the gauge symmetry, Wilson loops and
confinement are not. The `e^{-πβ}` scaling is specific to O(4); do not
transfer it to SU(2).

The O(3) paper's introduction cites Bałaban's work on lattice gauge theory and
Polyakov (1975); those are the primary sources to cite for the gauge-theory
side, rather than the O(3) paper itself.

## B. Mixing times and autocorrelations of heat-bath dynamics

Relevant to: [`burnin_autocorrelation.md`](burnin_autocorrelation.md),
[`monte_carlo.md`](monte_carlo.md), [`overrelaxation.md`](overrelaxation.md).

| Key | What it proves | Where to cite |
|---|---|---|
| `OAI:A-spectral-gap-throughout-the-high-temperature-Sherrington-Kirkpatrick-phase-September-24-2026` | For `β < 1`, the spectral gap of single-site heat-bath dynamics for the SK model stays bounded away from zero. | `burnin_autocorrelation.md`: "Why This Matters". An example of when single-site updates provably mix without critical slowing. |
| `OAI:Cutoff-throughout-the-high-temperature-Sherrington-Kirkpatrick-phase-September-24-2026` | Total-variation cutoff at `log n / (2λ(β))` for `β < 1`. | `burnin_autocorrelation.md`: burn-in length scales like `log(volume)` in the good regime. |
| `OAI:Critical-slowing-down-in-the-Sherrington-Kirkpatrick-model-September-24-2026` and `OAI:Critical-mixing-in-the-Sherrington-Kirkpatrick-model-September-25-2026` | At `β = 1`, mixing time `n^{2/3+o(1)}` (rate-one clocks), with a slow-mixing lower bound from typical equilibrium starts. | `overrelaxation.md`, "Why Overrelaxation Is Not Used Alone": a rigorous statement that local updates slow down polynomially at criticality. |
| `OAI:Functional-universality-of-critical-SK-autocorrelations-October-5-2026` and `OAI:Universality-of-critical-quench-autocorrelations-in-the-Sherrington-Kirkpatrick-model-October-5-2026` | Functional scaling limit of stationary and quench autocorrelations at `β = 1`, time scaled by `n^{2/3}`; the limit stays random. | `burnin_autocorrelation.md`, "Current Limitation": the integrated autocorrelation time is itself a random, sample-dependent quantity near criticality. |
| `OAI:Stretched-exponential-barriers-for-typical-SK-initial-states-September-24-2026` | For `β > 1`, stretched-exponential mixing barriers from typical equilibrium starts. | Optional. Only if the docs discuss frozen or metastable runs; the analogy to a confining gauge theory is weak. |

**Caveats.** The SK model is a mean-field spin glass with Gaussian couplings.
It has no lattice, no gauge symmetry and no locality, and its critical point
is a spin-glass transition, not a continuum limit. Use these papers only to
support general statements ("local-update Markov chains slow down near a
critical point, and the slowdown can be quantified in solvable models").
They give no numerical value of `z` or `τ_int` for lattice SU(2), and they
are not a substitute for the standard lattice references on autocorrelation
analysis (Sokal's lecture notes, Madras–Sokal, Wolff's Γ-method), which the
lab should cite directly.

## C. Spectral gaps of quantum Hamiltonians (methodology only)

Relevant to: [`spectroscopy.md`](spectroscopy.md), the finite-lattice
transfer-matrix framing.

| Key | What it proves | Where to cite |
|---|---|---|
| `OAI:The-periodic-spin-one-Haldane-gap-September-24-2026` | Uniform positive gap for the spin-1 Heisenberg chain on even periodic rings. | Optional, in a "related gap results" paragraph. |
| `OAI:Uniform-Stability-of-the-Spherical-Laughlin-Gap-October-5-2026` | Uniform gap for the Laughlin `V1` Hamiltonian at filling 1/3 on spheres, stable under weak perturbations. | Optional, same paragraph. |

**Caveat.** These are condensed-matter Hamiltonians with a gap proved
uniformly in system size. They show that "gap uniform in volume" is the
statement a rigorous mass-gap result needs; they say nothing about gauge
theories.

## D. Not recommended

- **Entropy "area law" papers** (family 265): the phrase collides with the
  Wilson area law. It is a statement about entanglement entropy and should
  not be cited next to string tension.
- **BFSS / matrix-model papers** (family 270): SU(N) matrix quantum mechanics
  is a dimensional reduction of Yang–Mills, but the preprint's result does not
  address the lattice theory studied here.

## A relevant check that is not a preprint

The same fork contains a note, `SU2_LATTICE_GAP.md`, which gives the exact
two-dimensional SU(2) string tension

    σ(β) = −ln( I₂(β) / I₁(β) )

and exactly solvable small systems. It is a fork-level computation, not part
of the OpenAI release, and it says so itself. Because this lab already has a
dimension-general SU(2) implementation, the 2D value is a candidate
validation target for the Wilson-loop and Creutz-ratio code in
`experiments/validation_suite.py`. I have not yet run that comparison.

## Suggested citation text

> Related rigorous work on two-dimensional asymptotically free models,
> including exact lattice mass asymptotics for the O(4) model and an
> O(3) continuum limit with a mass gap, is collected in the OpenAI Math
> Release preprints [OAI:…]. These do not concern four-dimensional gauge
> theory and do not bear on the present lattice study beyond providing
> context for what a rigorous gap statement looks like.
