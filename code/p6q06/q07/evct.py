#!/usr/bin/env python3
"""XVI-Q Q07 goal 1: sigma0-consistency of the evanescent counterterms that remove the ghost-row evanescent remainder.
CT_ev := - (evanescent part of every one-loop pole entering the rows: Gamma_{Kch}, Gamma_{cbar ch}, Gamma_{c* cc},
Gamma_{Kc}, Gamma_{cbar c}), i.e. the p22 prescription (local: each pole is a polynomial).  The rows are linear in the
loop inputs, so sigma0(Gamma_pole + CT_ev) on a row = the row evaluated on the 4D-PROJECTED inputs (every T set to 0 in
every loop input; tree vertices untouched, they are d-dimensional).  Test: that row vanishes IDENTICALLY (including any
T structure the tree vertices could generate).  If it does, the evanescent CTs are sigma0-consistent on these rows."""
from __future__ import annotations
import sys
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "p6q05"))
sys.path.insert(0, str(HERE.parent))
import rows as R  # noqa: E402
import rows6  # noqa: E402


def four_d(e):
    e = sp.sympify(e)
    return sp.expand(e.xreplace({x: 0 for x in e.free_symbols if x.name.startswith("T_")}))


def projected_parts(target):
    orig = (R.load, R.D13_own, R.D06_own, R.D13_K_ext)
    try:
        R.load = lambda g: four_d(orig[0](g))                       # noqa: E731
        R.D13_own = lambda *a: four_d(orig[1](*a))                  # noqa: E731
        R.D06_own = lambda *a: four_d(orig[2](*a))                  # noqa: E731
        R.D13_K_ext = lambda *a: four_d(orig[3](*a))                # noqa: E731
        return rows6.parts(target)
    finally:
        R.load, R.D13_own, R.D06_own, R.D13_K_ext = orig


def test(target):
    P = projected_parts(target)
    C = rows6.coeffs(target)
    Rr = sp.expand(sum(C[k] * P[k] for k in P))
    n_t = sum(1 for t in sp.Add.make_args(Rr) if any(x.name.startswith("T_") for x in t.free_symbols)) if Rr != 0 else 0
    # control: dropping the CT on one input (keep its evanescent pole) must re-open the row
    full = rows6.parts(target)
    ctrl = sp.expand(Rr + C["c1"] * (full["c1"] - P["c1"]))
    return {"row": target, "residual_with_evanescent_CTs": str(Rr), "zero": Rr == 0, "n_T_terms": n_t,
            "control_without_G10_evanescent_CT_nonzero": ctrl != 0}


if __name__ == "__main__":
    for t in ("ccK", "ccb"):
        r = test(t)
        print(t, r["zero"], r["n_T_terms"], r["control_without_G10_evanescent_CT_nonzero"], r["residual_with_evanescent_CTs"][:300])
