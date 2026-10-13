#!/usr/bin/env python3
"""XVI-Q Q09 (a): BRST-trivial families with DERIVED relative normalisation.
sigma(Y) = (S0, Y).  For Y = K.F(h,h): X_hhh = D-op(F) [dS0/dh . dY/dK], X_Khc = A-op(F) [dY/dh . s0 h].  The row
operators applied to X are again A-op and D-op; they commute, so sigma^2 = 0 on the rows forces
row(X) = A(D F) + lambda D(A F) = 0 with lambda = -1 (the antibracket's relative sign between the field and the
antifield derivative of a fermionic Y).  Same for Y = cbar.Phi(h,h) with the E-op (s cbar = sigma_F F(h)):
lambda_P = -1.  Y = c*.G(c,h): both pieces are dY/dc* . dS0/dc (same derivative type), lambda_G = +1 in the Q08
construction's normalisation.  Here: membership of every generator in N = ker(rows, Bose) at the derived lambda,
and the trace content of the failures."""
import json, pickle, sys  # noqa: E401
from pathlib import Path
import sympy as sp
from sympy.polys.matrices import DomainMatrix
HERE = Path(__file__).resolve().parent
Q = HERE.parents[4]
sys.path[:0] = [str(HERE.parents[3] / "p6q05"), str(HERE.parents[2]), str(HERE.parents[1]), str(HERE)]
from evfit import monomials  # noqa: E402
LOG = Q / "build" / "logs"
MS = {"K": monomials({"aK": 2, "gc1": 1, "a3": 2}, 3), "C": monomials({"gb2": 1, "gc1": 1, "a3": 2}, 4),
      "S": monomials({"z": 1, "gc1": 1, "gc2": 1}, 3), "H": monomials({"a1": 2, "a2": 2, "a3": 2}, 4)}
NAMES = [m for k in "KCSH" for m in MS[k]]
TRACE = {"H_aK_aK", "T_aK_aK", "H_a1_a1", "T_a1_a1", "H_a2_a2", "T_a2_a2", "H_a3_a3", "T_a3_a3"}
DERIVED = {"F": -1, "P": -1, "G": 1}


def has_trace(v):
    return any(c != 0 and {s.name for s in NAMES[i].free_symbols} & TRACE for i, c in enumerate(v))


def main():
    d = pickle.loads((LOG / "q08_triv.pkl").read_bytes())
    A = d["A"]
    res = {}
    vecs = []
    for name, pairs in d["fams"].items():
        lam = DERIVED[name]
        n_in, fails, fails_trace = 0, 0, 0
        for p, q in pairs:
            v = sp.Matrix(p) + lam * sp.Matrix(q)
            if all(sp.expand(x) == 0 for x in A * v):
                n_in += 1
                vecs.append(v)
            else:
                fails += 1
                fails_trace += has_trace(list(p) + []) or has_trace(list(q))
        res[name] = {"lambda_derived": str(lam), "n_generators": len(pairs), "n_in_N": n_in, "n_fail": fails,
                     "n_fail_with_trace_monomials": fails_trace}
    if vecs:
        res["dim_trivial_span_derived"] = DomainMatrix.from_Matrix(sp.Matrix.hstack(*vecs)).rank()
        Ag = d["Ag"]
        M = sp.Matrix.hstack(*vecs)
        res["dim_trivial_span_derived_inside_ghost_eq"] = res["dim_trivial_span_derived"] - \
            DomainMatrix.from_Matrix(Ag * M).rank()
    (LOG / "q09_fam.json").write_text(json.dumps(res, indent=1) + "\n")
    print(res)


if __name__ == "__main__":
    main()
