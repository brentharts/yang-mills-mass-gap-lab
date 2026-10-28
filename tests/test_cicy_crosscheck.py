"""Cross-check of the one-loop coefficients against CICY (brentharts/CICY).

The literal values below are the ones pyCICY.theories.running states
(MSSM_BETAS, the b_3 table for 3, 4 and 5 generations, and the Standard Model
b_3 = -7).  When CICY is importable (CICY_DIR or a sibling checkout) the same
numbers are also read from it.
"""
import math
import os
import sys
from pathlib import Path

import pytest

from ymlab import symbolic as S


def _cicy():
    root = os.environ.get("CICY_DIR") or str(Path(__file__).resolve().parents[2] / "CICY")
    if not (Path(root) / "pyCICY" / "theories" / "running.py").exists():
        return None
    sys.path.insert(0, root)
    try:
        from pyCICY.theories import running
    except Exception:
        return None
    return running


def test_susy_b3_matches_cicy_table_for_generations():
    b = S.one_loop_coefficient(susy=True)
    # n_f = 2 n_g flavour pairs for N = 3: b_3 = -9 + 2 n_g
    assert [b.numeric(3, 2 * ng) for ng in (3, 4, 5)] == [-3.0, -1.0, 1.0]


def test_standard_model_b3_is_minus_seven():
    assert S.one_loop_coefficient(susy=False).numeric(3, 6) == -7.0


def test_threshold_reproduces_four_generation_limit():
    # CICY: confines for n_g <= 4, not for n_g >= 5; n_f = 2 n_g, N = 3
    ineq = S.asymptotic_freedom_threshold(susy=True)
    limit = float(ineq.rhs.subs(S.N, 3))
    assert [2 * ng < limit for ng in (3, 4, 5)] == [True, True, False]


def test_lambda_exponent_is_the_lattice_spacing_exponent():
    assert S.lambda_exponent_identity()


def test_pure_glue_b0_is_the_coefficient_used_for_spacing():
    assert sp_eq(S.pure_glue_b0(), 11 * S.N / (48 * math.pi**0 * S.sp.pi**2))


def sp_eq(a, b):
    return S.sp.simplify(a - b) == 0


def test_live_cicy_agrees_when_available():
    running = _cicy()
    if running is None:
        pytest.skip("CICY not found")
    b = S.one_loop_coefficient(susy=True)
    for ng in (3, 4, 5):
        assert float(running.beta_coefficients(ng, 1)[2]) == b.numeric(3, 2 * ng)
    assert running.confines(4) and not running.confines(5)
