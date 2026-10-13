#!/usr/bin/env python3
"""XVI-Q Q11 diagnostic: for every Q08 family generator, which block of the linear system (row ccK, row ccb,
row hhc, Bose symmetry of X_hhh) the two pieces violate."""
import json, pickle, sys  # noqa: E401
from pathlib import Path
import sympy as sp
HERE = Path(__file__).resolve().parent
Q = HERE.parents[6]
sys.path[:0] = [str(HERE.parents[5] / "p6q05"), str(HERE.parents[4]), str(HERE.parents[3]), str(HERE.parents[2])]
import rows as R  # noqa: E402
import evfit2 as E2  # noqa: E402
from evfit import ansatz  # noqa: E402
LOG = Q / "build" / "logs"


def blocks():
    XK, xk = ansatz({"aK": 2, "gc1": 1, "a3": 2}, 3, "xk")
    XC, xc = ansatz({"gb2": 1, "gc1": 1, "a3": 2}, 4, "xc")
    XS, xs = ansatz({"z": 1, "gc1": 1, "gc2": 1}, 3, "xs")
    XH, xh = ansatz({"a1": 2, "a2": 2, "a3": 2}, 4, "xh")
    sig = sp.sympify(json.loads((LOG / "q07_hhc_tree.json").read_text())["sigma_F"])
    xall = xk + xc + xs + xh
    X = set(xall)
    out = {}
    rows = {"ccK": E2.ghost_row_total("ccK", XK, XS) - E2.ghost_row_total("ccK", sp.Integer(0), sp.Integer(0)),
            "ccb": E2.ghost_row_total("ccb", XC, XS) - E2.ghost_row_total("ccb", sp.Integer(0), sp.Integer(0)),
            "hhc": E2.hhc_ansatz_total(XH, XK, XC, sig) - E2.hhc_ansatz_total(sp.Integer(0), sp.Integer(0), sp.Integer(0), sig)}
    for k, r in rows.items():
        r = sp.expand(r)
        gens = sorted([s for s in r.free_symbols if s not in X], key=str)
        out[k] = sp.linear_eq_to_matrix(sp.Poly(r, *gens).coeffs(), xall)[0]
    beq = []
    for mp in ({"p1": R.P2, "p2": R.P1, "a1": "a2", "a2": "a1"}, {"p1": R.P3, "p2": R.P2, "a1": "a3", "a3": "a1"}):
        d = sp.expand(XH - R.remap(XH, mp))
        gens = sorted([s for s in d.free_symbols if s not in X], key=str)
        beq += sp.Poly(d, *gens).coeffs()
    out["bose"] = sp.linear_eq_to_matrix(beq, xall)[0]
    return out


def main():
    B = blocks()
    pickle.dump(B, open(LOG / "q11_blocks.pkl", "wb"))
    d = pickle.loads((LOG / "q08_triv.pkl").read_bytes())
    res = {}
    for name, pairs in d["fams"].items():
        cnt = {}
        for p, q in pairs:
            P, Qv = sp.Matrix(p), sp.Matrix(q)
            if not any(P) and not any(Qv):
                continue
            key = []
            for b, M in B.items():
                ap = any(sp.expand(x) != 0 for x in M * P)
                aq = any(sp.expand(x) != 0 for x in M * Qv)
                key.append(f"{b}:{'P' if ap else ''}{'Q' if aq else ''}")
            k = " ".join(key)
            cnt[k] = cnt.get(k, 0) + 1
        res[name] = cnt
    (LOG / "q11_diag.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
