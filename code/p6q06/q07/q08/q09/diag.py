#!/usr/bin/env python3
"""XVI-Q Q09 diagnostic: evaluate the rows on one BRST-trivial family vector, row by row."""
import json, pickle, sys  # noqa: E401
from pathlib import Path
import sympy as sp
HERE = Path(__file__).resolve().parent
Q = HERE.parents[4]
sys.path[:0] = [str(HERE.parents[3] / "p6q05"), str(HERE.parents[2]), str(HERE.parents[1]), str(HERE.parent), str(HERE)]
import evfit2 as E2  # noqa: E402
from evfit import monomials  # noqa: E402
LOG = Q / "build" / "logs"
I = sp.I
MS = {"K": monomials({"aK": 2, "gc1": 1, "a3": 2}, 3), "C": monomials({"gb2": 1, "gc1": 1, "a3": 2}, 4),
      "S": monomials({"z": 1, "gc1": 1, "gc2": 1}, 3), "H": monomials({"a1": 2, "a2": 2, "a3": 2}, 4)}


def split(v):
    out, o = {}, 0
    for k in "KCSH":
        n = len(MS[k])
        out[k] = sp.expand(I * sum(c * m for c, m in zip(v[o:o + n], MS[k])))
        o += n
    return out


def rows(v):
    X = split(v)
    sig = sp.sympify(json.loads((LOG / "q07_hhc_tree.json").read_text())["sigma_F"])
    base = sp.sympify(json.loads((LOG / "q07_hhc_row.json").read_text())["residual"])
    r = {"ccK": E2.ghost_row_total("ccK", X["K"], X["S"]) - E2.ghost_row_total("ccK", 0, 0),
         "ccb": E2.ghost_row_total("ccb", X["C"], X["S"]) - E2.ghost_row_total("ccb", 0, 0),
         "hhc": E2.hhc_ansatz_total(X["H"], X["K"], X["C"], sig) - base}
    return {k: sp.expand(v_) for k, v_ in r.items()}, X


if __name__ == "__main__":
    fam, lam, k = sys.argv[1], sp.sympify(sys.argv[2]), int(sys.argv[3])
    d = pickle.loads((LOG / "q08_triv.pkl").read_bytes())
    p, q = d["fams"][fam][k]
    v = [sp.sympify(a) + lam * sp.sympify(b) for a, b in zip(p, q)]
    r, X = rows(v)
    print({kk: len(sp.Add.make_args(x)) if x != 0 else 0 for kk, x in r.items()})
    print({kk: str(x)[:300] for kk, x in X.items() if x != 0})
    print({kk: str(x)[:300] for kk, x in r.items() if x != 0})
