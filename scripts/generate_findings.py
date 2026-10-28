"""
Generate the findings markdown, the LaTeX paper, and (optionally) the PDF.

    python scripts/generate_findings.py [--author NAME] [--no-pdf]

All numbers come from results/runs, results/findings and ymlab.symbolic; the
figures from scripts/generate_figures.py are regenerated first.
"""
import argparse
import numpy as np
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from ymlab import figures, findings, symbolic  # noqa: E402

FIGDIR = ROOT / "results" / "figures" / "findings"
PAPER = ROOT / "paper"
CITE = ("Elsbernd, \\emph{Yang--Mills Mass Gap Laboratory} "
        "\\url{https://github.com/Elsbernd09/yang-mills-mass-gap-lab}")


def sched(a):
    return a["schedule"]


def agg_rows(d, tex):
    out = []
    for a in d["aggregation"]:
        cells = [f"{a['beta']:g}", sched(a).replace("+", "+") if not tex else sched(a),
                 f"{a['mean_gain']:.2f}", f"{a['mean_cost_factor']:.3f}",
                 f"{a['mean_of_ratios']:.3f}", f"{a['ratio_of_means']:.3f}",
                 f"{a['min']:.2f}--{a['max']:.2f}" if tex else f"{a['min']:.2f}–{a['max']:.2f}"]
        out.append(cells)
    return out


def md_table(head, rows):
    s = "| " + " | ".join(head) + " |\n|" + "---|" * len(head) + "\n"
    return s + "\n".join("| " + " | ".join(r) + " |" for r in rows) + "\n"


def tex_table(head, rows, align):
    s = f"\\begin{{tabular}}{{{align}}}\n\\toprule\n" + " & ".join(head) + " \\\\\n\\midrule\n"
    s += "\n".join(" & ".join(r) + " \\\\" for r in rows)
    return s + "\n\\bottomrule\n\\end{tabular}\n"


HEAD = ["$\\beta$", "schedule", "mean gain", "mean cost factor",
        "mean of ratios", "ratio of means", "pair range"]


def repro_text(d, tex):
    if not d["reproduction"]:
        return "No local reproduction file was present when this document was generated."
    r = d["reproduction"][0]
    esc = (lambda x: str(x).replace("_", "\\_")) if tex else str
    le = {k: esc(v) for k, v in r["local_env"].items()}
    pe = {k: esc(v) for k, v in r["published_env"].items()}
    t = (f"On {le['platform_system']} {le['platform_machine']} (Python {le['python_version']}, "
         f"NumPy {le['numpy_version']}) the published seed-2026 conditions were re-run; the published "
         f"runs were made on {pe['platform_system']} {pe['platform_machine']} (Python "
         f"{pe['python_version']}, NumPy {pe['numpy_version']}). The pure-Metropolis chains reproduce the "
         f"published $\\tau_{{int}}$ to within {r['metropolis_tau_max_diff']:.1e}; "
         f"{r['hybrid_identical']} of {r['hybrid_total']} hybrid chains reproduce it exactly, with local-to-published "
         f"$\\tau_{{int}}$ ratios between {r['hybrid_min_tau_ratio']:.2f} and {r['hybrid_max_tau_ratio']:.2f}. "
         f"Median runtime here was {r['runtime_slowdown_median']:.1f}$\\times$ the published runtime. "
         "The cause of the hybrid divergence was not isolated (platform, library version, and floating-point "
         "ordering are all confounded).")
    return t if tex else t.replace("$\\tau_{int}$", "τ_int").replace("$\\times$", "×")



def on_text(c):
    """(O(N) line, SU(2) line, note) for the comparison section, or None."""
    if c is None:
        return None
    on = "; ".join(
        "O(%d): measured/predicted local slope %s over xi = %.1f to %.1f (smallest L/xi %.1f)" % (
            n, ", ".join("%.2f" % x for x in r["ratio"]), r["xi_range"][0], r["xi_range"][1],
            r["L_over_xi_min"])
        for n, r in sorted(c["on"].items()))
    su = []
    worst = []
    for size, r in sorted(c["su"].items()):
        if "ratio" in r:
            su.append("%d^4 lattice: measured/two-loop slope %s at beta = %s" % (
                size, ", ".join("%.2f" % x for x in r["ratio"]),
                ", ".join("%.1f" % x for x in r["mid"])))
            worst.append(float(np.max(r["ratio"])))
    dev = "; ".join("%d^4: largest |plaquette - (1 - 3/4beta)| at beta >= 2.8 is %.3f" % kv
                    for kv in sorted(c["plaquette_weak_dev"].items()))
    sizes = sorted(k for k, r in c["su"].items() if "ratio" in r)
    grow = ""
    if len(sizes) >= 2:
        a, b = c["su"][sizes[0]]["ratio"][-1], c["su"][sizes[-1]]["ratio"][-1]
        grow = (" (at the largest beta the ratio is %.2f on %d^4 and %.2f on %d^4, "
                "so the larger lattice %s)" % (a, sizes[0], b, sizes[-1],
                "moved it toward one" if b > a else "did not move it toward one"))
    note = ("The SU(2) lattices are far from asymptotic scaling: the slope of ln sqrt(chi) is "
            "%.0f%% of the two-loop value at best, and at most %.0f%% over the last three betas of every lattice%s, a possible "
            "reason is strong-coupling and finite-size effects in chi(2,2) on a lattice this small, which was not tested here. The plaquette, by contrast, is already close to its weak-coupling leading term "
            "at beta of 2.8 and above (%s). The O(N) models reach %.2f to %.2f of the asymptotic "
            "slope at their largest beta (correlation lengths up to %.0f lattice spacings); the SU(2) scans are not "
            "run at lattice sizes where a correlation length could be measured at all."
            % (100 * max(worst), 100 * max(max(r['ratio'][-3:]) for r in c['su'].values() if 'ratio' in r), grow, dev, min(r['ratio'][-1] for r in c['on'].values()),
               max(r['ratio'][-1] for r in c['on'].values()),
               max(r['xi_range'][1] for r in c['on'].values()))) if worst else "No SU(2) scan files were present."
    return on, "; ".join(su) or "No SU(2) scan files were present.", note


def on_md():
    c = findings.on_comparison()
    if c is None:
        return ""
    on_line, su_line, su_note = on_text(c)
    f2 = {r["model"]: r["delta_beta"] for r in c["factor_two"]}
    return f"""## 6. Asymptotic freedom: 2D O(N) against 4D SU(2)

The physics tools of the OpenAI-math fork ([brentharts/openai_math_yang_mills](https://github.com/brentharts/openai_math_yang_mills), commit `{c['fork_commit'][:9]}`) include a Wolff cluster Monte Carlo of the two-dimensional O(N) model; its output is stored in `results/findings/openai_fork_on_model_mc.json`. The leading logarithmic slope of the mass (O(N)) or lattice spacing (SU(N)) against beta is derived with sympy (`ymlab.symbolic.on_mass_log_slope`, `su_n_log_spacing_slope`). For SU(N) it is -12 pi^2 / (11 N^2) plus a correction 51 / (121 beta) that does not depend on N.

**O(N), measured against the form.** {on_line}. The measured slope approaches the asymptotic form from below as beta grows, which is the usual approach and needed correlation lengths of order ten lattice spacings.

![O(N) local slopes](../results/figures/findings/on_local_slopes.png)

**How far in beta a scale change takes.** Halving the scale takes about {f2['O(3), beta=1.5']:.2f} in beta for O(3), {f2['O(4), beta=2.0']:.2f} for O(4), {f2['SU(2), beta=2.5']:.2f} for SU(2) and {f2['SU(3), beta=6.0']:.2f} for SU(3), at the couplings named.

**SU(2) in four dimensions** (`scripts/su2_scaling_scan.py`, Metropolis, one chain per beta). {su_line}

![SU(2) scan](../results/figures/findings/su2_scan.png)
![SU(2) slopes](../results/figures/findings/su2_slopes.png)

{su_note}

*Limits.* The O(N) models are a toy for Yang-Mills, with different dimension, observables and exponents. The SU(2) lattices are tiny, short, and use a single chain per beta, with errors blocked inside that chain. The O(N) length is a second-moment estimate and not the transfer-matrix gap. Nothing here is a statement about the continuum mass gap.

"""


def on_tex():
    c = findings.on_comparison()
    if c is None:
        return ""
    on_line, su_line, su_note = on_text(c)
    f2 = {r["model"]: r["delta_beta"] for r in c["factor_two"]}
    esc = lambda t: (t.replace("^4", "$^4$").replace("%", "\\%").replace("_", "\\_")
                      .replace("sqrt(chi)", "$\\sqrt{\\chi}$").replace("beta", "$\\beta$")
                      .replace("chi(2,2)", "$\\chi(2,2)$"))
    return (
        "\\section{Asymptotic freedom: two-dimensional O(N) against four-dimensional SU(2)}\n"
        "The physics tools of the OpenAI-math fork "
        "(\\url{https://github.com/brentharts/openai_math_yang_mills}, commit \\texttt{%s}) include a Wolff "
        "cluster Monte Carlo of the two-dimensional O(N) model, whose output is stored with the laboratory. "
        "The leading logarithmic slope of the mass (O(N)) or lattice spacing (SU(N)) against $\\beta$ is "
        "derived with sympy; for SU(N) it is $-12\\pi^2/(11N^2)$ plus a correction $51/(121\\beta)$ "
        "independent of $N$.\n\n"
        "\\paragraph{O(N).} %s. The measured slope approaches the asymptotic form from below as $\\beta$ "
        "grows, which needed correlation lengths of order ten lattice spacings.\n"
        "\\begin{figure}[ht]\\centering\\includegraphics[width=0.85\\linewidth]{on_local_slopes}"
        "\\caption{Measured local slope of $\\ln m$ against $\\beta$ and the asymptotic form.}\\end{figure}\n"
        "Halving the scale takes $\\Delta\\beta\\approx %.2f$ for O(3), %.2f for O(4), %.2f for SU(2) and %.2f "
        "for SU(3) at the couplings named in the text file.\n\n"
        "\\paragraph{SU(2).} %s.\n"
        "\\begin{figure}[ht]\\centering\\includegraphics[width=0.85\\linewidth]{su2_scan}"
        "\\caption{Plaquette and Creutz ratio from the SU(2) scans; dashed and dotted lines are the leading "
        "strong- and weak-coupling terms.}\\end{figure}\n"
        "\\begin{figure}[ht]\\centering\\includegraphics[width=0.6\\linewidth]{su2_slopes}"
        "\\caption{Measured slope of $\\ln\\sqrt{\\chi}$ against the two-loop SU(2) slope.}\\end{figure}\n"
        "%s\n\n"
        "These comparisons are descriptive: the O(N) models are a toy for Yang--Mills, the SU(2) lattices are "
        "tiny with one chain per $\\beta$, and nothing here concerns the continuum mass gap.\n\n" % (
            c["fork_commit"][:9], esc(on_line), f2["O(3), beta=1.5"], f2["O(4), beta=2.0"],
            f2["SU(2), beta=2.5"], f2["SU(3), beta=6.0"], esc(su_line), esc(su_note)))


def cicy_md():
    from ymlab import symbolic as S
    b = S.one_loop_coefficient(susy=True)
    row = ", ".join("n_g = %d: %+.0f" % (ng, b.numeric(3, 2 * ng)) for ng in (3, 4, 5))
    return f"""## 7. One-loop cross-check with CICY

[CICY](https://github.com/brentharts/CICY) (`pyCICY.theories.running`) computes the one-loop coefficients b_i of a supersymmetric spectrum from an index and uses b_3 for the QCD scale, Lambda = M exp(-2 pi / (|b_3| alpha)). With the convention d(1/alpha)/d ln mu = b / (2 pi), sympy gives b = -(11/3)N + (2/3) n_f for Dirac fundamentals and b = -3N + n_f for chiral flavour pairs (`ymlab.symbolic.one_loop_coefficient`). For N = 3 and n_f = 2 n_g flavour pairs this gives SU(3) b_3 = {row}, the values CICY states, and b_3 = -7 for the Standard Model. Asymptotic freedom then needs n_f < 3N, i.e. at most four generations, which is CICY's `confines` statement. Sympy also shows that CICY's exponent 2 pi / (|b| alpha) is the same expression as the exponent 1 / (2 b_0 g^2) of the lattice spacing used in section 6, with b_0 = -b / (16 pi^2); at n_f = 0 that is 11N / (48 pi^2). The check is of the formulas against each other (`tests/test_cicy_crosscheck.py`); it adds no evidence about the lattice data.

"""


def cicy_tex():
    from ymlab import symbolic as S
    b = S.one_loop_coefficient(susy=True)
    row = ", ".join("$n_g=%d$: $%+.0f$" % (ng, b.numeric(3, 2 * ng)) for ng in (3, 4, 5))
    return (
        "\\section{One-loop cross-check with CICY}\n"
        "The CICY package (\\url{https://github.com/brentharts/CICY}, \\texttt{pyCICY.theories.running}) "
        "computes one-loop coefficients $b_i$ from an index and uses $b_3$ in the QCD scale "
        "$\\Lambda = M\\exp(-2\\pi/(|b_3|\\alpha))$. With $d(1/\\alpha)/d\\ln\\mu = b/2\\pi$, sympy gives "
        "$b=-\\tfrac{11}{3}N+\\tfrac{2}{3}n_f$ (Dirac fundamentals) and $b=-3N+n_f$ (supersymmetric pairs). "
        "For $N=3$, $n_f=2n_g$ this reproduces SU(3) $b_3$: " + row + ", the values CICY states, and "
        "$b_3=-7$ for the Standard Model; asymptotic freedom needs $n_f<3N$, at most four generations. "
        "The exponent $2\\pi/(|b|\\alpha)$ equals the lattice-spacing exponent $1/(2b_0g^2)$ "
        "with $b_0=-b/16\\pi^2$, which is $11N/48\\pi^2$ at $n_f=0$. This compares formulas with formulas; "
        "it adds no evidence about the lattice data.\n\n")


def build_md(d):
    det = d["determinism"]
    m = d["manifest"]["configuration"]
    s = f"""# Findings: efficiency decomposition, aggregation, and reproducibility

*Generated by `scripts/generate_findings.py` from `{d['campaign']}`, `results/findings/` and `ymlab.symbolic`. Do not edit by hand.*

Builds on Michael Elsbernd's Yang–Mills Mass Gap Laboratory (`ymlab`) and its paired overrelaxation campaign (β ∈ {m['betas']}, seeds {m['seeds']}, lattice {m['shape']}, {m['measurement_updates']} measured updates). This is a **finite-lattice computational-methodology** study; nothing here is a statement about the continuum Yang–Mills mass gap.

## 1. Efficiency ratio decomposes into gain × cost

With N_eff = N/(2τ) and η = N_eff/t, the paired ESS/s ratio of a hybrid schedule to Metropolis is

    R = η_OR / η_M = (τ_M / τ_OR) · (t_M / t_OR)

The number of measurements cancels. R > 1 requires the decorrelation gain τ_M/τ_OR to exceed the cost penalty t_OR/t_M. Sympy derives this (`ymlab.symbolic.efficiency_ratio`) and the test-suite checks it against the recorded ESS/s ratios of every pair to 1e-9.

![gain vs cost](../results/figures/findings/gain_vs_cost.png)

**Finding.** Across all {d['n_pairs']} paired (β, seed, schedule) conditions the plaquette ratio was between {d['ratio_min']:.2f} and {d['ratio_max']:.2f}; {d['n_ratio_above_one']} exceeded 1. The gain exceeded 1 in {d['gains_above_one']} of {d['n_pairs']} pairs (max {d['gain_max']:.2f}) but the cost penalty was between {d['penalty_min']:.1f} and {d['penalty_max']:.1f}. In this configuration the extra overrelaxation sweeps were never repaid in plaquette ESS/s.

![ESS by schedule](../results/figures/findings/ess_by_schedule.png)

## 2. Aggregation rule changes the reported number

"Mean of paired ratios" and "ratio of mean ESS/s" are different statistics (the mean of products is not the product of means). Both are tabulated; gain and cost factor are per-pair means.

{md_table([h.replace('$\\beta$', 'β') for h in HEAD], agg_rows(d, False))}
The two rules differ by up to {100 * d['max_aggregation_gap']:.0f}% of the ratio of means in this table. With three seeds per cell, per-pair ratios are widely dispersed and neither summary should be read as more than descriptive.

![aggregation](../results/figures/findings/aggregation_sensitivity.png)

## 3. What is deterministic and what is not

Two published campaigns with the same configuration hash were compared: maximum |Δτ_int| = {det['tau_max_abs_diff']:.1e}, maximum relative ESS/s difference = {100 * det['ess_max_rel_diff']:.1f}% (same machine). τ_int is deterministic; ESS/s inherits wall-clock noise.

![determinism](../results/figures/findings/determinism.png)

## 4. Reproduction on another platform

{repro_text(d, False)}

The stored data is in `results/findings/`; rerun with `python scripts/reproduction_check.py`.

## 5. Symbolic models (sympy → LaTeX)

`src/ymlab/symbolic.py` models the lattice coupling β = 2N/g², the SU(2) plaquette term 1 − a₀, trace conjugation invariance, N_eff, the efficiency ratio and its break-even, the AR(1) τ_int = (1+φ)/(2(1−φ)), the cosh effective-mass identity, and the Creutz-ratio cancellation of perimeter and constant terms in an area law. Each has a numeric check against the numerical code in `tests/test_symbolic.py`. Illustrations: `ar1_tau`, `area_law_creutz`, `effective_mass_plateau` in `results/figures/findings/`.

{on_md()}{cicy_md()}## Limits

Three seeds, one lattice shape (5×4×4), SU(2), one observable emphasised (plaquette). The ESS/s ranking is specific to this implementation's Python/NumPy cost structure. Hybrid-chain divergence across platforms is observed, not explained.

## Reproduce

    python scripts/generate_figures.py
    python scripts/reproduction_check.py     # ~5 min
    python scripts/generate_findings.py      # writes this file and paper/
"""
    return s


def build_tex(d, author):
    det = d["determinism"]
    m = d["manifest"]["configuration"]
    fig = lambda n, cap, w="0.8": (
        f"\\begin{{figure}}[ht]\\centering\\includegraphics[width={w}\\linewidth]{{{n}}}"
        f"\\caption{{{cap}}}\\end{{figure}}\n")
    body = f"""\\documentclass[11pt]{{article}}
\\usepackage[margin=1in]{{geometry}}
\\usepackage{{amsmath,amssymb,booktabs,graphicx,hyperref}}
\\title{{Efficiency decomposition, aggregation sensitivity and reproducibility\\\\
in a finite-lattice overrelaxation study}}
\\author{{{author}}}
\\date{{\\today}}
\\begin{{document}}
\\maketitle
\\begin{{abstract}}
We analyse the paired overrelaxation-efficiency campaign of the Yang--Mills Mass Gap
Laboratory of Elsbernd. The per-pair ESS-per-second ratio factors exactly into a
decorrelation gain and a cost factor; in all {d['n_pairs']} paired conditions the cost
penalty ({d['penalty_min']:.1f}--{d['penalty_max']:.1f}$\\times$) outweighed the gain
(at most {d['gain_max']:.2f}$\\times$), so no ratio exceeded one (maximum {d['ratio_max']:.2f}).
We show that two natural aggregation rules disagree, that $\\tau_{{\\mathrm{{int}}}}$ is deterministic while
ESS/s is not, and report a cross-platform re-run in which hybrid chains did not reproduce.
This is a finite-lattice methodology study and makes no statement about the continuum mass gap.
\\end{{abstract}}

\\section{{Setting}}
The study uses the lattice and campaign of {CITE}: SU(2) Wilson action on a
${'\\times'.join(str(x) for x in m['shape'])}$ lattice, $\\beta\\in\\{{{', '.join(f'{b:g}' for b in m['betas'])}\\}}$,
seeds {', '.join(str(x) for x in m['seeds'])}, {m['measurement_updates']} measured updates per chain, and the
schedules M, M+1OR, M+2OR, M+4OR (one Metropolis sweep followed by 0, 1, 2, 4 microcanonical overrelaxation sweeps).
The campaign and its manifests are Elsbernd's; the analysis below only post-processes
the recorded results.

\\section{{Symbolic models}}
The following relations were derived or verified with \\texttt{{sympy}} (\\texttt{{ymlab.symbolic}}) and are
numerically cross-checked against the laboratory code in the test-suite.
{symbolic.latex_table()}

\\section{{Efficiency decomposition}}
{fig('gain_vs_cost', 'Per-pair decorrelation gain against cost penalty (log axes). Points below the dashed line would have $R>1$; none do.', '0.62')}
Across all {d['n_pairs']} pairs the plaquette ratio lay in [{d['ratio_min']:.2f}, {d['ratio_max']:.2f}];
the gain exceeded one in {d['gains_above_one']} pairs (maximum {d['gain_max']:.2f}) while the cost penalty lay in
[{d['penalty_min']:.1f}, {d['penalty_max']:.1f}]. For this implementation the additional sweeps were not repaid.
{fig('ess_by_schedule', 'Plaquette ESS/s by schedule and $\\\\beta$; dots are individual seeds.', '0.62')}

\\section{{Aggregation sensitivity}}
The mean of paired ratios and the ratio of mean ESS/s are different statistics.
\\begin{{table}}[ht]\\centering\\small
{tex_table(HEAD, agg_rows(d, True), 'rrrrrrr')}\\caption{{Aggregates per $(\\beta,\\text{{schedule}})$, three seeds each.}}
\\end{{table}}
They differ by up to {100 * d['max_aggregation_gap']:.0f}\\% of the ratio of means.
{fig('aggregation_sensitivity', 'Per-seed ratios with the two aggregates; dashed line is break-even.', '0.95')}

\\section{{Determinism and reproduction}}
Two published campaigns with identical configuration hash agree in $\\tau_{{\\mathrm{{int}}}}$ to
{det['tau_max_abs_diff']:.1e} (maximum) and in ESS/s to within {100 * det['ess_max_rel_diff']:.1f}\\%.
{fig('determinism', 'Run-to-run agreement on one machine.', '0.85')}
{repro_text(d, True)}

{on_tex()}{cicy_tex()}\\section{{Limits}}
Three seeds, one lattice, SU(2), plaquette emphasised; the ranking depends on the
Python/NumPy cost structure. The cross-platform divergence is observed, not explained.

\\paragraph{{Acknowledgement.}} All simulation code, the campaign design and the recorded data are the work
of Michael Elsbernd; this document is an analysis layer on top of that work.
\\end{{document}}
"""
    return body


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--author", default="ymlab analysis contributors")
    ap.add_argument("--no-pdf", action="store_true")
    a = ap.parse_args()
    figures.generate_static(FIGDIR)
    c = findings.on_comparison()
    if c is not None:
        from ymlab import scaling_comparison as SC
        figures.fig_on_local_slopes(c["on"], FIGDIR)
        figures.fig_factor_two(c["factor_two"], FIGDIR)
        if c["scans"]:
            figures.fig_su2_scan(c["scans"], FIGDIR)
            figures.fig_su2_slopes(c["scans"], FIGDIR)
    d = findings.collect()
    (ROOT / "docs" / "rosettamath_integration_findings.md").write_text(build_md(d))
    PAPER.mkdir(exist_ok=True)
    for p in FIGDIR.glob("*.pdf"):
        shutil.copy(p, PAPER / p.name)
    tex = PAPER / "ymlab_findings.tex"
    tex.write_text(build_tex(d, a.author))
    if not a.no_pdf:
        for _ in range(2):
            subprocess.run(["pdflatex", "-interaction=nonstopmode", tex.name],
                           cwd=PAPER, check=True, stdout=subprocess.DEVNULL)
        print(PAPER / "ymlab_findings.pdf")


if __name__ == "__main__":
    main()
