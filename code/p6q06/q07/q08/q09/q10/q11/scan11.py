#!/usr/bin/env python3
"""XVI-Q Q11 diagnostic (not used for any claim of triviality): per F generator, the set of (lam, mu) with
rows(D F + mu E phi + lam (A F + mu A phi)) = 0, to locate the normalisation that the derivation must reproduce."""
import json, pickle, sys  # noqa: E401
from pathlib import Path
import sympy as sp
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import sigma11 as S  # noqa: E402
L, M = sp.symbols("lam mu")


def main():
    S.SIG = sp.sympify(json.loads((S.LOG / "q07_hhc_tree.json").read_text())["sigma_F"])
    A = pickle.loads((S.LOG / "q08_triv.pkl").read_bytes())["A"]
    sols = {}
    for m in S.monomials({"aK": 2, "a2": 2, "a3": 2}, 2):
        f = sp.expand(m + S.R.remap(m, S.SW23))
        if f == 0:
            continue
        phi = S.cbar_partner(f)
        dF = sp.Matrix(S.vec("H", S.D_op(f))[0])
        aF = sp.Matrix(S.vec("K", S.A_on_K(f))[0])
        eP = sp.Matrix(S.vec("H", S.E_op(phi))[0]) if phi != 0 else sp.zeros(S.NDIM, 1)
        aP = sp.Matrix(S.vec("C", S.A_on_C(phi))[0]) if phi != 0 else sp.zeros(S.NDIM, 1)
        r = (A * dF + M * (A * eP) + L * (A * aF) + L * M * (A * aP)).applyfunc(sp.expand)
        eqs = [e for e in r if e != 0]
        if not eqs:
            k = "any"
        else:
            s = sp.solve(eqs, [L, M], dict=True)
            k = str(sorted([(str(x.get(L, "free")), str(x.get(M, "free"))) for x in s])) if s else "none"
        sols[k] = sols.get(k, 0) + 1
    (S.LOG / "q11_scan.json").write_text(json.dumps(sols, indent=1) + "\n")
    print(json.dumps(sols, indent=1))


if __name__ == "__main__":
    main()
