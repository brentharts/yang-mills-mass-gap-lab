"""
4D SU(2) beta scan for the scaling comparison with the 2D O(N) models.

For each (size, beta): thermalise with Metropolis, then measure the plaquette
and the planar Wilson loops W(1,1), W(2,1), W(2,2) (averaged over all six
planes) every `stride` sweeps.  Stores means with a blocked-jackknife error and
the Creutz ratio chi(2,2) = -ln[W22 W11 / (W21 W12)] (with W(1,2)=W(2,1) by
symmetry of the plane average) in results/findings/su2_scan_<L>.json.

    python scripts/su2_scaling_scan.py --size 4 --therm 200 --meas 600
"""
import argparse
import json
import platform
import sys
from itertools import combinations
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from ymlab.lattice import Lattice  # noqa: E402
from ymlab.monte_carlo import metropolis_sweep  # noqa: E402
from ymlab.observables import average_rectangular_wilson_loop as W  # noqa: E402
from ymlab.plaquette import average_plaquette  # noqa: E402

PLANES = list(combinations(range(4), 2))


def planar(lat, w, h):
    vals = []
    for mu, nu in PLANES:
        vals.append(W(lat, mu, nu, w, h))
        if w != h:
            vals.append(W(lat, mu, nu, h, w))
    return float(np.mean(vals))


def chi_series(w11, w21, w22):
    return -np.log(w22 * w11 / (w21 * w21))


def one(args):
    size, beta, therm, meas, stride, seed, eps = args
    lat = Lattice((size,) * 4, cold_start=True, seed=seed)
    for _ in range(therm):
        metropolis_sweep(lat, beta, eps)
    rows = []
    acc = []
    for i in range(meas):
        acc.append(metropolis_sweep(lat, beta, eps))
        if (i + 1) % stride == 0:
            rows.append([average_plaquette(lat), planar(lat, 1, 1),
                         planar(lat, 2, 1), planar(lat, 2, 2)])
    a = np.array(rows)
    nb = 10
    blocks = np.array_split(a, nb)
    bm = np.array([b.mean(axis=0) for b in blocks])
    chis = []
    for k in range(nb):  # jackknife over blocks
        m = np.delete(bm, k, axis=0).mean(axis=0)
        chis.append(chi_series(m[1], m[2], m[3]))
    mean = a.mean(axis=0)
    chi = float(chi_series(mean[1], mean[2], mean[3]))
    jk = float(np.sqrt((nb - 1) / nb * np.sum((np.array(chis) - np.mean(chis)) ** 2)))
    return {"beta": beta, "size": size, "plaquette": float(mean[0]),
            "w11": float(mean[1]), "w21": float(mean[2]), "w22": float(mean[3]),
            "chi22": chi, "chi22_err": jk, "acceptance": float(np.mean(acc)),
            "therm": therm, "meas": meas, "stride": stride, "seed": seed}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--size", type=int, default=4)
    ap.add_argument("--therm", type=int, default=200)
    ap.add_argument("--meas", type=int, default=600)
    ap.add_argument("--stride", type=int, default=2)
    ap.add_argument("--eps", type=float, default=0.3)
    ap.add_argument("--betas", type=float, nargs="+", default=[2.0, 2.2, 2.4, 2.6, 2.8, 3.0])
    ap.add_argument("--procs", type=int, default=1)
    a = ap.parse_args()
    jobs = [(a.size, b, a.therm, a.meas, a.stride, 2026 + i, a.eps) for i, b in enumerate(a.betas)]
    with Pool(a.procs) as p:
        res = p.map(one, jobs)
    out = {"environment": {"system": platform.system(), "machine": platform.machine(),
                           "python": platform.python_version(), "numpy": np.__version__},
           "points": res}
    dest = ROOT / "results" / "findings"
    dest.mkdir(parents=True, exist_ok=True)
    path = dest / f"su2_scan_{a.size}.json"
    path.write_text(json.dumps(out, indent=2))
    print(path)


if __name__ == "__main__":
    main()
