#!/usr/bin/env python3
"""XVI-Q Q11 goal 2 (conditional): decomposition of N_g (rows + Bose + ghost equation) into
 (i) the span of the Ktilde-based sigma-images that pass (derived signs), (ii) blind coordinates, (iii) the rest."""
import json, pickle  # noqa: E401
from pathlib import Path
import sympy as sp
from sympy.polys.matrices import DomainMatrix
HERE = Path(__file__).resolve().parent
LOG = HERE.parents[6] / "build" / "logs"


def rk(M):
    return DomainMatrix.from_Matrix(M).rank() if M.shape[1] else 0


def main():
    d = pickle.loads((LOG / "q08_triv.pkl").read_bytes())
    A, Ag = d["A"], d["Ag"]
    n = A.shape[1]
    keep = pickle.loads((LOG / "q11_sigma_vecs.pkl").read_bytes())
    V = keep[("F", "derived_minus1")] + keep[("G", "derived_minus1")]
    M = sp.Matrix.hstack(*[sp.Matrix(v) for v in V])
    zc = [j for j in range(n) if all(A[i, j] == 0 for i in range(A.shape[0]))]
    Bm = sp.Matrix.hstack(*[sp.Matrix([1 if i == j else 0 for i in range(n)]) for j in zc])
    Bg = Bm.shape[1] - rk(Ag * Bm)                       # blind directions inside the ghost-eq subspace
    full = A.col_join(Ag)
    dimNg = n - rk(full)
    r_triv = rk(M)
    r_tb = rk(M.row_join(Bm))
    # blind directions satisfying the ghost equation: kernel of Ag restricted to blind coords
    ker = (Ag * Bm).nullspace()
    Bk = sp.Matrix.hstack(*[Bm * k for k in ker]) if ker else sp.zeros(n, 0)
    r_tbg = rk(M.row_join(Bk))
    res = {"dim_N_g": dimNg, "dim_sigma_span_passing": r_triv, "dim_blind_in_N_g": Bg,
           "dim_sigma_plus_blindg": r_tbg, "dim_sigma_cap_blindg": r_triv + Bg - r_tbg,
           "dim_remaining_N_g_mod_(sigma+blind)": dimNg - r_tbg,
           "n_passing_vectors": {"F": len(keep[("F", "derived_minus1")]), "G": len(keep[("G", "derived_minus1")])},
           "dim_sigma_plus_all_blind": r_tb}
    (LOG / "q11_class.json").write_text(json.dumps(res, indent=1) + "\n")
    print(res)


if __name__ == "__main__":
    main()
