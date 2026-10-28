import math

import numpy as np
import sympy as sp

from ymlab import scaling_comparison as SC
from ymlab import symbolic as S


def _synthetic_on(n, betas):
    xi = [1.0 / (b ** (1 / (n - 2)) * math.exp(-2 * math.pi * b / (n - 2))) for b in betas]
    pts = [{"beta": b, "xi": x, "L_over_xi": 9.0} for b, x in zip(betas, xi)]
    return {"data": {"results": {str(n): pts}}}


def test_on_slope_matches_numerical_derivative():
    m = S.on_mass_log_slope()
    for n in (3, 4, 5):
        b, h = 1.7, 1e-6
        f = lambda x: math.log(x ** (1 / (n - 2)) * math.exp(-2 * math.pi * x / (n - 2)))
        assert abs((f(b + h) - f(b - h)) / (2 * h) - m.numeric(n, b)) < 1e-6


def test_su_n_leading_coefficient_and_nindependent_correction():
    m = S.su_n_log_spacing_slope()
    lead = lambda n: -12 * math.pi**2 / (11 * n**2)
    for n in (2, 3):
        assert abs(m.numeric(n, 1e6) - lead(n)) < 1e-5
    # the one-over-beta correction is the same number for every N
    c2 = m.numeric(2, 3.0) - lead(2)
    c3 = m.numeric(3, 3.0) - lead(3)
    assert abs(c2 - c3) < 1e-12 and abs(c2 - 51 / (121 * 3.0)) < 1e-12


def test_su2_leading_slope_is_3pi2_over_11():
    assert abs(-12 * math.pi**2 / (11 * 4) + 3 * math.pi**2 / 11) < 1e-12


def test_local_slopes_recover_exact_exponential():
    b = np.array([1.0, 1.5, 2.0])
    mid, s = SC.local_slopes(b, -3.0 * b)
    assert np.allclose(s, -3.0) and np.allclose(mid, [1.25, 1.75])


def test_on_rows_on_exact_scaling_data_gives_ratio_near_one():
    rows = SC.on_rows(_synthetic_on(4, [1.6, 1.8, 2.0, 2.2]))
    r = rows[4]["ratio"]
    assert np.all(np.abs(r - 1) < 0.02)


def test_factor_two_table_ordering():
    rows = {r["model"].split(",")[0]: r["delta_beta"] for r in SC.factor_two_table()}
    assert rows["O(3)"] < rows["O(4)"] < rows["SU(2)"] < rows["SU(3)"]


def test_su2_rows_handles_nonpositive_chi():
    scan = {"points": [
        {"beta": 2.0, "chi22": 0.5, "chi22_err": 0.05, "plaquette": 0.5},
        {"beta": 2.2, "chi22": -0.1, "chi22_err": 0.05, "plaquette": 0.6},
        {"beta": 2.4, "chi22": 0.2, "chi22_err": 0.05, "plaquette": 0.65},
    ]}
    r = SC.su2_rows(scan)
    assert list(r["valid"]) == [True, False, True]
    assert len(r["measured"]) == 1


def test_plaquette_limits():
    assert SC.strong_coupling_plaquette(2.0) == 0.5
    assert abs(SC.weak_coupling_plaquette(3.0) - 0.75) < 1e-15
