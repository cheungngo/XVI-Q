#!/usr/bin/env python3
"""XVI-Q Q11 diagnostic: trace content of F generators vs. membership of the derived sigma-image in N."""
import json, pickle, sys  # noqa: E401
from pathlib import Path
import sympy as sp
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import sigma11 as S  # noqa: E402
TR = {"H_aK_aK", "H_a2_a2", "H_a3_a3", "T_aK_aK", "T_a2_a2", "T_a3_a3"}


def main():
    S.SIG = sp.sympify(json.loads((S.LOG / "q07_hhc_tree.json").read_text())["sigma_F"])
    A = pickle.loads((S.LOG / "q08_triv.pkl").read_bytes())["A"]
    tab = {}
    for m in S.monomials({"aK": 2, "a2": 2, "a3": 2}, 2):
        f = sp.expand(m + S.R.remap(m, S.SW23))
        if f == 0:
            continue
        names = {s.name for s in m.free_symbols}
        tr = "trace" if names & TR else "no_trace"
        ev = "evanescent" if any(n.startswith("T_") for n in names) else "4d"
        phi = S.cbar_partner(f)
        v = S.add(S.vec("H", S.D_op(f))[0], S.vec("H", S.E_op(phi))[0] if phi != 0 else [0] * S.NDIM,
                  S.scale(-1, S.add(S.vec("K", S.A_on_K(f))[0], S.vec("C", S.A_on_C(phi))[0] if phi != 0 else [0] * S.NDIM)))
        v0 = S.add(S.vec("H", S.D_op(f))[0], S.scale(-1, S.vec("K", S.A_on_K(f))[0]))
        ok = all(sp.expand(x) == 0 for x in A * sp.Matrix(v))
        ok0 = all(sp.expand(x) == 0 for x in A * sp.Matrix(v0))
        k = f"{tr}/{ev}: Ktilde={'in' if ok else 'out'} K-only={'in' if ok0 else 'out'}"
        tab[k] = tab.get(k, 0) + 1
    (S.LOG / "q11_trace.json").write_text(json.dumps(tab, indent=1, sort_keys=True) + "\n")
    print(json.dumps(tab, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
