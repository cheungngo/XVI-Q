#!/usr/bin/env python3
"""XVI-Q Q10: anatomy of the Q08/Q09 BRST-family vectors (reads the frozen q08_triv.pkl; no new loop data).
For every generator (P = first piece, Q = second piece) of F, G, P:
  * blind  : the generator's vector at the derived lambda is supported only on coordinates that the rows+Bose
             system A does not see at all (zero columns of A), i.e. membership in N is vacuous;
  * kpure  : support only in the K sector (X_Khc);
  * in_N   : A(P + lam Q) = 0 at the derived lambda;
  * lam_sol: the set of lambda making A(P + lam Q) = 0 (None = no lambda; 'any' = both pieces in N).
Also: dims of the blind subspace of N, of the derived trivial span modulo blind directions, and the >=36 bound."""
import json, pickle, sys  # noqa: E401
from pathlib import Path
import sympy as sp
from sympy.polys.matrices import DomainMatrix
HERE = Path(__file__).resolve().parent
Q = HERE.parents[5]
sys.path[:0] = [str(HERE.parents[4] / "p6q05"), str(HERE.parents[3]), str(HERE.parents[2]), str(HERE.parents[1]), str(HERE)]
from evfit import monomials  # noqa: E402
LOG = Q / "build" / "logs"
MS = {"K": monomials({"aK": 2, "gc1": 1, "a3": 2}, 3), "C": monomials({"gb2": 1, "gc1": 1, "a3": 2}, 4),
      "S": monomials({"z": 1, "gc1": 1, "gc2": 1}, 3), "H": monomials({"a1": 2, "a2": 2, "a3": 2}, 4)}
SECT = [k for k in "KCSH" for _ in MS[k]]
DERIVED = {"F": -1, "P": -1, "G": 1}
L = sp.Symbol("lam")


def rk(M):
    return DomainMatrix.from_Matrix(M).rank() if M.shape[1] else 0


def main():
    d = pickle.loads((LOG / "q08_triv.pkl").read_bytes())
    A, Ag = d["A"], d["Ag"]
    n = A.shape[1]
    zero_cols = {j for j in range(n) if all(A[i, j] == 0 for i in range(A.shape[0]))}
    res = {"n_unknowns": n, "n_blind_coordinates": len(zero_cols),
           "blind_coordinates_by_sector": {k: sum(1 for j in zero_cols if SECT[j] == k) for k in "KCSH"}}
    passing, passing_nonblind = [], []
    for name, pairs in d["fams"].items():
        lam = DERIVED[name]
        st = {"n": len(pairs), "in_N": 0, "in_N_blind": 0, "in_N_nonblind": 0, "fail": 0, "fail_kpure": 0,
              "zero_vector": 0, "lam_any": 0, "lam_unique": {}, "lam_none": 0}
        for p, q in pairs:
            P, Qv = sp.Matrix(p), sp.Matrix(q)
            v = P + lam * Qv
            supp = {j for j in range(n) if v[j] != 0}
            if not supp:
                st["zero_vector"] += 1
                continue
            blind = supp <= zero_cols
            kpure = all(SECT[j] == "K" for j in supp)
            r = (A * v).applyfunc(sp.expand)
            if all(x == 0 for x in r):
                st["in_N"] += 1
                passing.append(v)
                if blind:
                    st["in_N_blind"] += 1
                else:
                    st["in_N_nonblind"] += 1
                    passing_nonblind.append(v)
            else:
                st["fail"] += 1
                st["fail_kpure"] += kpure
            ap, aq = (A * P).applyfunc(sp.expand), (A * Qv).applyfunc(sp.expand)
            if all(x == 0 for x in ap) and all(x == 0 for x in aq):
                st["lam_any"] += 1
                continue
            eqs = [sp.expand(a + L * b) for a, b in zip(ap, aq) if sp.expand(a + L * b) != 0]
            sol = sp.solve(eqs, L, dict=True)
            if sol:
                key = str(sol[0][L])
                st["lam_unique"][key] = st["lam_unique"].get(key, 0) + 1
            else:
                st["lam_none"] += 1
        res[name] = st
    Ad = A
    res["dim_N"] = n - rk(Ad)
    res["dim_N_with_ghost_eq"] = n - rk(A.col_join(Ag))
    # blind subspace of N: spanned by unit vectors on zero columns
    res["dim_blind_subspace"] = len(zero_cols)
    Bmat = sp.Matrix.hstack(*[sp.Matrix([1 if i == j else 0 for i in range(n)]) for j in sorted(zero_cols)]) \
        if zero_cols else sp.zeros(n, 0)
    res["dim_blind_inside_ghost_eq"] = Bmat.shape[1] - rk(Ag * Bmat) if Bmat.shape[1] else 0
    if passing:
        M = sp.Matrix.hstack(*passing)
        res["dim_trivial_span_derived"] = rk(M)
        res["dim_trivial_span_derived_mod_blind"] = rk(M.row_join(Bmat)) - Bmat.shape[1]
        res["dim_trivial_span_derived_inside_ghost_eq"] = rk(M) - rk(Ag * M)
        MB = M.row_join(Bmat)
        res["dim_trivial_plus_blind"] = rk(MB)
    (LOG / "q10_blind.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
