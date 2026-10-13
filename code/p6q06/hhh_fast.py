#!/usr/bin/env python3
"""XVI-Q Q06 P.6 hhh graphs with the fast ring path (fastpole.py):
  G01  V3^3 triangle (one orientation, symmetry 1; momenta as p6q05/graphs.G01)
  G02  V3 (one external h) -- two internal gravitons -- V4 (two external h), symmetry 1/2, summed over the 3 choices
usage: hhh_fast.py G01 | G02_1 | G02_2 | G02_3 | G02 (sum) -> build/logs/q06_<id>.json"""
from __future__ import annotations
import json, sys, time  # noqa: E401
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
Q = HERE.parents[1]
sys.path.insert(0, str(HERE.parent / "p6q05"))
sys.path.insert(0, str(HERE))
import gravity as G  # noqa: E402
import graphs as GR  # noqa: E402
import uvpole as U  # noqa: E402
import fastpole as FP  # noqa: E402
import v4  # noqa: E402

p1, p2, l, p3 = GR.p1, GR.p2, GR.l, GR.p3
EXT = {"1": ("a1", p1), "2": ("a2", p2), "3": ("a3", p3)}
W = 4
LOG = Q / "build" / "logs"


def _names(probes):
    return ["a1", "a2", "a3", "p1", "p2", "l"] + [f"{x}{y}" for x in probes for y in "UW"]


def G01():
    F = FP.Fast(_names(["hAB", "hAC", "hBA", "hBC", "hCB", "hCA"]))
    A = F.trunc(F.to(GR._v3("a1", p1, "hAB", -(l + p1), "hAC", l)), W)
    B = F.trunc(F.to(GR._v3("a2", p2, "hBA", l + p1, "hBC", -(l + p1 + p2))), W)
    C = F.trunc(F.to(GR._v3("a3", p3, "hCB", l + p1 + p2, "hCA", -l)), W)
    n = F.contract(A, B, [("h", "hAB", "hBA")], W)
    n = F.contract(n, C, [("h", "hBC", "hCB"), ("h", "hAC", "hCA")], W)
    return GR.r4(F.loop_pole(n, [U.V(), p1, p1 + p2], W, 1), "G01")


def G02(ch):
    a, k = EXT[ch]
    o = [EXT[c] for c in "123" if c != ch]
    F = FP.Fast(_names(["hA1", "hA2", "hB1", "hB2"]))
    A = F.trunc(F.to(GR._v3(a, k, "hA1", -l, "hA2", l - k)), W)
    B = F.trunc(F.to(v4.V4_grav((GR.ext(o[0][0]), o[0][1]), (GR.ext(o[1][0]), o[1][1]),
                                (G.probe_mat("hB1"), l), (G.probe_mat("hB2"), k - l))), W)
    n = F.contract(A, B, [("h", "hA1", "hB1"), ("h", "hA2", "hB2")], W)
    return GR.r4(F.loop_pole(n, [U.V(), -k], W, sp.Rational(1, 2)), f"G02 ch{ch}")


def run(gid):
    t0 = time.time()
    if gid == "G02":
        r = sp.expand(sum(sp.sympify(json.loads((LOG / f"q06_G02_{c}.json").read_text())["pole_over_U_d4_residue"])
                          for c in "123"))
    elif gid == "G01":
        r = G01()
    else:
        r = G02(gid[-1])
    out = {"graph": "P6-" + gid, "pole_over_U_d4_residue": str(r), "n_terms": len(sp.Add.make_args(r)),
           "zero": r == 0, "seconds": round(time.time() - t0, 1)}
    (LOG / f"q06_{gid}.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


if __name__ == "__main__":
    o = run(sys.argv[1])
    print(sys.argv[1], o["n_terms"], o["zero"], o["seconds"])
