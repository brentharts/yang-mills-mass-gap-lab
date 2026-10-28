"""
Asymptotic-freedom scaling: 2D O(N) Monte Carlo against 4D SU(2) lattice data.

The O(N) numbers come from the physics tools of the OpenAI-math fork
(brentharts/openai_math_yang_mills, tools/on_model_mc.py), stored here under
results/findings/ with the fork commit they came from.  The SU(2) numbers come
from scripts/su2_scaling_scan.py.  Both are compared with the slopes derived
symbolically in ymlab.symbolic.

Descriptive only.  The O(N) models are a toy for Yang-Mills, the SU(2) lattices
here are tiny, and neither says anything about the continuum mass gap.
"""

from __future__ import annotations

import json
import math
import os
import subprocess
from pathlib import Path

import numpy as np

from . import symbolic as S

ROOT = Path(__file__).resolve().parents[2]
FINDINGS = ROOT / "results" / "findings"
FORK_FILE = FINDINGS / "openai_fork_on_model_mc.json"


def locate_fork() -> Path | None:
    for c in (os.environ.get("OPENAI_MATH_DIR"), ROOT.parent / "openai_math_yang_mills",
              Path.home() / "openai_math_yang_mills"):
        if c and (Path(c) / "tools" / "on_model_mc.py").exists():
            return Path(c)
    return None


def import_fork(dest: Path = FORK_FILE) -> Path:
    """Copy the fork's build/data/on_model_mc.json next to the provenance."""
    fork = locate_fork()
    if fork is None:
        raise FileNotFoundError("set OPENAI_MATH_DIR to the fork checkout")
    src = fork / "build" / "data" / "on_model_mc.json"
    payload = json.loads(src.read_text())
    commit = subprocess.run(["git", "-C", str(fork), "rev-parse", "HEAD"],
                            capture_output=True, text=True).stdout.strip()
    remote = subprocess.run(["git", "-C", str(fork), "remote", "get-url", "origin"],
                            capture_output=True, text=True).stdout.strip()
    out = {"source": "tools/on_model_mc.py via make models", "fork_remote": remote,
           "fork_commit": commit, "data": payload}
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=2))
    return dest


def load_on(path: Path = FORK_FILE) -> dict | None:
    path = Path(path)
    return json.loads(path.read_text()) if path.exists() else None


def local_slopes(beta: np.ndarray, log_value: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Midpoint betas and slopes d(log_value)/d(beta) between adjacent points."""
    beta = np.asarray(beta, float)
    log_value = np.asarray(log_value, float)
    return 0.5 * (beta[1:] + beta[:-1]), np.diff(log_value) / np.diff(beta)


def on_rows(on: dict) -> dict[int, dict]:
    """Per N: betas, xi, local slopes of ln(1/xi), predicted slope, ratio."""
    model = S.on_mass_log_slope()
    rows = {}
    for n, pts in on["data"]["results"].items():
        n = int(n)
        b = np.array([p["beta"] for p in pts])
        xi = np.array([p["xi"] for p in pts])
        mid, meas = local_slopes(b, -np.log(xi))  # ln m = -ln xi
        pred = np.array([model.numeric(n, m) for m in mid])
        rows[n] = {"beta": b, "xi": xi, "mid": mid, "measured": meas,
                   "predicted": pred, "ratio": meas / pred,
                   "xi_range": (float(xi.min()), float(xi.max())),
                   "L_over_xi_min": min(p["L_over_xi"] for p in pts)}
    return rows


def su2_rows(scan: dict) -> dict:
    """Local slopes of ln sqrt(chi22) against beta from an SU(2) scan file."""
    pts = sorted(scan["points"], key=lambda p: p["beta"])
    b = np.array([p["beta"] for p in pts])
    chi = np.array([p["chi22"] for p in pts])
    err = np.array([p["chi22_err"] for p in pts])
    plaq = np.array([p["plaquette"] for p in pts])
    ok = chi > 0
    model = S.su_n_log_spacing_slope()
    out = {"beta": b, "chi": chi, "chi_err": err, "plaquette": plaq, "valid": ok}
    if ok.sum() >= 2:
        mid, slope = local_slopes(b[ok], 0.5 * np.log(chi[ok]))
        # error from the two endpoints' jackknife errors (independent chains)
        rel = 0.5 * err[ok] / chi[ok]
        serr = np.sqrt(rel[1:] ** 2 + rel[:-1] ** 2) / np.diff(b[ok])
        pred = np.array([model.numeric(2, m) for m in mid])
        out.update({"mid": mid, "measured": slope, "slope_err": serr,
                    "predicted": pred, "ratio": slope / pred})
    return out


def strong_coupling_plaquette(beta: float) -> float:
    """Leading strong-coupling term <(1/2) Tr U_p> = beta/4 for SU(2)."""
    return beta / 4.0


def weak_coupling_plaquette(beta: float) -> float:
    """Leading weak-coupling term <(1/2) Tr U_p> = 1 - 3/(4 beta) for SU(2)."""
    return 1.0 - 3.0 / (4.0 * beta)


def factor_two_table() -> list[dict]:
    """Change in beta that halves the scale, at the betas each model is run."""
    on = S.on_mass_log_slope()
    su = S.su_n_log_spacing_slope()
    rows = []
    for label, slope in (
        ("O(3), beta=1.5", abs(on.numeric(3, 1.5))),
        ("O(4), beta=2.0", abs(on.numeric(4, 2.0))),
        ("SU(2), beta=2.5", abs(su.numeric(2, 2.5))),
        ("SU(3), beta=6.0", abs(su.numeric(3, 6.0))),
    ):
        rows.append({"model": label, "slope": slope,
                     "delta_beta": S.beta_for_factor_two(slope)})
    return rows
