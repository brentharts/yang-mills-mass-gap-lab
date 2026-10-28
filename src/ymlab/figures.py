"""Publication figures for the integration findings (Agg backend, PDF+PNG)."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from . import efficiency_analysis as EA  # noqa: E402
from .spectroscopy import arccosh_effective_mass, periodic_cosh_correlator  # noqa: E402

# colour-blind-safe (Okabe-Ito subset)
C = {"M": "#000000", 1: "#0072B2", 2: "#E69F00", 4: "#D55E00", "ref": "#7f7f7f"}
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.25,
                     "axes.spines.top": False, "axes.spines.right": False})


def _save(fig, outdir: Path, name: str) -> list[Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    paths = []
    for ext in ("pdf", "png"):
        p = outdir / f"{name}.{ext}"
        fig.savefig(p, dpi=200, bbox_inches="tight")
        paths.append(p)
    plt.close(fig)
    return paths


def fig_ess_by_schedule(records, outdir):
    betas = sorted({r["beta"] for r in records})
    sweeps = sorted({int(r["overrelaxation_sweeps"]) for r in records})
    fig, ax = plt.subplots(figsize=(5.2, 3.2))
    w = 0.2
    for j, s in enumerate(sweeps):
        for i, b in enumerate(betas):
            v = [r["plaquette_ess_per_second"] for r in records
                 if r["beta"] == b and int(r["overrelaxation_sweeps"]) == s]
            x = i + (j - 1.5) * w
            ax.bar(x, np.mean(v), w, color=C["M" if s == 0 else s], alpha=0.55,
                   label=("M" if s == 0 else f"M+{s}OR") if i == 0 else None)
            ax.plot([x] * len(v), v, "o", ms=3, color=C["M" if s == 0 else s])
    ax.set_xticks(range(len(betas)), [f"{b:g}" for b in betas])
    ax.set_xlabel(r"$\beta$"); ax.set_ylabel("plaquette ESS / s")
    ax.legend(frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.15))
    return _save(fig, outdir, "ess_by_schedule")


def fig_gain_vs_cost(rows, outdir):
    fig, ax = plt.subplots(figsize=(4.8, 3.6))
    g = np.geomspace(0.5, 10, 50)
    ax.plot(g, g, color=C["ref"], ls="--", lw=1)
    ax.text(3.0, 3.6, r"break-even $R=1$", color=C["ref"], rotation=30, fontsize=8)
    for s in (1, 2, 4):
        sel = [r for r in rows if r["sweeps"] == s]
        ax.scatter([r["gain"] for r in sel], [r["penalty"] for r in sel],
                   color=C[s], s=22, label=f"M+{s}OR")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"decorrelation gain $\tau_M/\tau_{OR}$")
    ax.set_ylabel(r"cost penalty $t_{OR}/t_M$")
    ax.legend(frameon=False)
    return _save(fig, outdir, "gain_vs_cost")


def fig_aggregation(rows, records, outdir):
    agg = EA.aggregate(rows); rom = EA.ratio_of_means(records)
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.8), sharey=True)
    betas = sorted({a["beta"] for a in agg})
    for ax, b in zip(axes, betas):
        for k, s in enumerate((1, 2, 4)):
            sel = [r["ratio"] for r in rows if r["beta"] == b and r["sweeps"] == s]
            ax.plot([k] * len(sel), sel, "o", ms=3, color=C[s], alpha=0.6)
            a = next(x for x in agg if x["beta"] == b and x["sweeps"] == s)
            ax.plot(k - 0.15, a["mean_of_ratios"], "_", ms=14, mew=2, color="k")
            ax.plot(k + 0.15, rom[(b, s)], "D", ms=4, mfc="none", color="k")
        ax.axhline(1, color=C["ref"], ls="--", lw=1)
        ax.set_xticks(range(3), ["+1", "+2", "+4"]); ax.set_title(rf"$\beta={b:g}$")
    axes[0].set_ylabel("ESS/s ratio vs Metropolis")
    axes[1].plot([], [], "_", color="k", mew=2, label="mean of paired ratios")
    axes[1].plot([], [], "D", mfc="none", color="k", label="ratio of means")
    axes[1].legend(frameon=False, fontsize=7, loc="upper right")
    return _save(fig, outdir, "aggregation_sensitivity")


def fig_determinism(rec_a, rec_b, outdir):
    key = lambda r: (r["beta"], r["seed"], r["overrelaxation_sweeps"])
    b = {key(r): r for r in rec_b}
    pairs = [(r, b[key(r)]) for r in rec_a]
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.9))
    for ax, k, ttl in ((axes[0], "plaquette_tau_int", r"$\tau_{int}$ (deterministic)"),
                       (axes[1], "plaquette_ess_per_second", "ESS/s (timing-dependent, same machine)")):
        x = [p[0][k] for p in pairs]; y = [p[1][k] for p in pairs]
        ax.loglog(x, y, "o", ms=3, color=C[1])
        lo, hi = min(x + y), max(x + y)
        ax.plot([lo, hi], [lo, hi], color=C["ref"], ls="--", lw=1)
        ax.set_title(ttl); ax.set_xlabel("campaign run A"); ax.set_ylabel("campaign run B")
    return _save(fig, outdir, "determinism")


def fig_area_law(outdir, sigma=0.3, perimeter=0.2, const=0.1):
    r = np.arange(2, 7)
    W = lambda a, b: np.exp(-sigma * a * b - perimeter * (a + b) - const)
    chi = [-np.log(W(i, i) * W(i - 1, i - 1) / (W(i, i - 1) * W(i - 1, i))) for i in r]
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.8))
    for i in range(1, 6):
        axes[0].semilogy(range(1, 7), [W(i, t) for t in range(1, 7)], "o-", ms=3, label=f"R={i}")
    axes[0].set_xlabel("T"); axes[0].set_ylabel("W(R,T)"); axes[0].legend(frameon=False, fontsize=7)
    axes[1].plot(r, chi, "o", color=C[1]); axes[1].axhline(sigma, color=C["ref"], ls="--")
    axes[1].set_ylim(0, 0.5); axes[1].set_xlabel("R=T"); axes[1].set_ylabel(r"$\chi(R,R)$")
    axes[1].text(2, sigma + 0.02, r"$\sigma$", color=C["ref"])
    return _save(fig, outdir, "area_law_creutz")


def fig_effective_mass(outdir, mass=0.7, T=16):
    t = np.arange(0, T + 1)
    c = periodic_cosh_correlator(t, 1.0, mass, T)
    rng = np.random.default_rng(1)
    noisy = c * (1 + 0.01 * rng.standard_normal(len(c)))
    fig, ax = plt.subplots(figsize=(4.6, 3.0))
    ax.plot(t, arccosh_effective_mass(c), "-", color=C["ref"], label="exact")
    ax.plot(t, arccosh_effective_mass(noisy), "o", ms=3, color=C[1], label="1% noise")
    ax.axhline(mass, color="k", ls=":", lw=1)
    ax.set_xlabel("t"); ax.set_ylabel(r"$m_{\rm eff}(t)$"); ax.legend(frameon=False)
    return _save(fig, outdir, "effective_mass_plateau")


def fig_ar1_tau(outdir):
    ph = np.linspace(-0.5, 0.95, 200)
    fig, ax = plt.subplots(figsize=(4.4, 3.0))
    ax.semilogy(ph, (1 + ph) / (2 * (1 - ph)), color=C[1])
    ax.axhline(0.5, color=C["ref"], ls="--", lw=1)
    ax.set_xlabel(r"AR(1) coefficient $\varphi$"); ax.set_ylabel(r"$\tau_{int}$")
    return _save(fig, outdir, "ar1_tau")


def fig_reproduction(published, local, outdir):
    """published/local: dict (beta)->(plaquette_ess_per_second M, M+1OR)."""
    betas = sorted(local)
    fig, ax = plt.subplots(figsize=(4.8, 3.0))
    x = np.arange(len(betas))
    ax.bar(x - 0.2, [published[b] for b in betas], 0.4, color=C[ "ref"], label="published run")
    ax.bar(x + 0.2, [local[b] for b in betas], 0.4, color=C[1], label="reproduced here")
    ax.axhline(1, color="k", lw=0.8, ls="--")
    ax.set_xticks(x, [f"{b:g}" for b in betas]); ax.set_xlabel(r"$\beta$")
    ax.set_ylabel("M+1OR / M  ESS-per-second"); ax.legend(frameon=False)
    return _save(fig, outdir, "reproduction")


def generate_static(outdir: Path, campaign: Path | None = None) -> dict[str, list[Path]]:
    """Generate every figure that needs only the published campaign data."""
    outdir = Path(outdir)
    dirs = EA.campaign_dirs()
    chosen = Path(campaign) if campaign else dirs[-1]
    rec = EA.load_records(chosen)
    rows = EA.paired_decomposition(rec)
    out = {
        "ess_by_schedule": fig_ess_by_schedule(rec, outdir),
        "gain_vs_cost": fig_gain_vs_cost(rows, outdir),
        "aggregation_sensitivity": fig_aggregation(rows, rec, outdir),
        "area_law_creutz": fig_area_law(outdir),
        "effective_mass_plateau": fig_effective_mass(outdir),
        "ar1_tau": fig_ar1_tau(outdir),
    }
    if len(dirs) >= 2:
        out["determinism"] = fig_determinism(
            EA.load_records(dirs[0]), EA.load_records(dirs[1]), outdir)
    return out


# --- asymptotic-freedom comparison -----------------------------------------
def fig_on_local_slopes(on_rows, outdir):
    """Measured local slope of ln m against beta, over the one-loop form."""
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.9))
    for ax, (n, r) in zip(axes, sorted(on_rows.items())):
        ax.plot(r["mid"], r["predicted"], "-", color=C["ref"], label="asymptotic form")
        ax.plot(r["mid"], r["measured"], "o", color=C[1], label="Monte Carlo")
        ax.set_title(f"O({n}), 2D")
        ax.set_xlabel(r"$\beta$"); ax.set_ylabel(r"$d\ln m/d\beta$")
        ax.annotate("ratio " + ", ".join(f"{x:.2f}" for x in r["ratio"]),
                    (0.04, 0.06), xycoords="axes fraction", fontsize=7)
    axes[0].legend(frameon=False, fontsize=7, loc="upper right")
    return _save(fig, outdir, "on_local_slopes")


def fig_su2_scan(scans, outdir):
    """Plaquette and Creutz ratio against beta for each scanned size."""
    from . import scaling_comparison as SC
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.9))
    markers = {4: "o", 6: "s"}
    bs = np.linspace(1.9, 3.1, 50)
    axes[0].plot(bs, [SC.strong_coupling_plaquette(b) for b in bs], "--", color=C["ref"], lw=1,
                 label=r"$\beta/4$")
    axes[0].plot(bs, [SC.weak_coupling_plaquette(b) for b in bs], ":", color=C["ref"], lw=1,
                 label=r"$1-3/4\beta$")
    for size, scan in sorted(scans.items()):
        r = SC.su2_rows(scan)
        axes[0].plot(r["beta"], r["plaquette"], markers.get(size, "^"), ms=4, color=C[1 if size == 4 else 4],
                     label=f"{size}$^4$")
        v = r["valid"]
        axes[1].errorbar(r["beta"][v], r["chi"][v], r["chi_err"][v], fmt=markers.get(size, "^"), ms=4,
                         color=C[1 if size == 4 else 4], label=f"{size}$^4$", capsize=2)
    axes[0].set_xlabel(r"$\beta$"); axes[0].set_ylabel("plaquette"); axes[0].legend(frameon=False, fontsize=7)
    axes[1].set_xlabel(r"$\beta$"); axes[1].set_ylabel(r"Creutz ratio $\chi(2,2)$")
    return _save(fig, outdir, "su2_scan")


def fig_su2_slopes(scans, outdir):
    from . import scaling_comparison as SC
    fig, ax = plt.subplots(figsize=(4.6, 3.0))
    for size, scan in sorted(scans.items()):
        r = SC.su2_rows(scan)
        if "mid" not in r:
            continue
        ax.errorbar(r["mid"], r["measured"], r["slope_err"], fmt="os"[size == 6], ms=4,
                    color=C[1 if size == 4 else 4], capsize=2, label=f"{size}$^4$ measured")
    from . import symbolic as S
    m = S.su_n_log_spacing_slope()
    bs = np.linspace(2.0, 3.0, 20)
    ax.plot(bs, [m.numeric(2, b) for b in bs], "-", color=C["ref"], label="two-loop SU(2)")
    ax.set_xlabel(r"$\beta$"); ax.set_ylabel(r"$d\ln\sqrt{\chi}/d\beta$"); ax.legend(frameon=False, fontsize=7)
    return _save(fig, outdir, "su2_slopes")


def fig_factor_two(rows, outdir):
    fig, ax = plt.subplots(figsize=(4.4, 2.8))
    ax.barh([r["model"] for r in rows], [r["delta_beta"] for r in rows], color=C[1])
    ax.set_xlabel(r"change in $\beta$ that halves the scale")
    ax.invert_yaxis()
    return _save(fig, outdir, "factor_two")
