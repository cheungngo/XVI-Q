"""XVI-Q Q01 O.2(a): INDEPENDENT reconstruction of the full 10-component ell_4 and exact comparison with the XVI-P p23
export (3360 terms).  Route (independent of XVI-P's quartic.py): the Gamma-Gamma Lagrangian
    L_GG = sqrt(-g) g^{mn} (Gamma^a_{mb} Gamma^b_{na} - Gamma^a_{mn} Gamma^b_{ab}),   g = eta + eps H,
with H_{mk}, dH_{mk,r} independent jet symbols, g^{-1} by the Neumann series and sqrt(-g) by exp(tr log / 2), both
truncated exactly at eps^4.  ell_4 := s * [eps^4] L_GG with the overall sign/normalisation s fixed by matching XVI-O's
conformal probe ONLY (one number), then the full polynomial is compared term by term."""
from __future__ import annotations
import json
import sys
import time
from pathlib import Path

import sympy as sp

QROOT = Path(__file__).resolve().parents[2]
ETA = [1, -1, -1, -1]
EPS = sp.Symbol("eps")


def Hs(m, k):
    a, b = sorted((m, k))
    return sp.Symbol(f"H{a}{b}")


def dHs(m, k, r):
    a, b = sorted((m, k))
    return sp.Symbol(f"dH{a}{b}_{r}")


def build():
    """exact eps-grading: [eps^4] L_GG = sum over (H-degree split) with two Gammas (each one dH); computed in the
    sparse polynomial ring QQ[H, dH] (exact)."""
    from sympy.polys.rings import ring
    from sympy.polys.domains import QQ
    names = [f"H{a}{b}" for a in range(4) for b in range(a, 4)] + \
            [f"dH{a}{b}_{r}" for a in range(4) for b in range(a, 4) for r in range(4)]
    R, *gens = ring(",".join(names), QQ)
    gd = dict(zip(names, gens))
    H = [[gd[f"H{min(m, k)}{max(m, k)}"] for k in range(4)] for m in range(4)]
    dH = [[[gd[f"dH{min(m, k)}{max(m, k)}_{r}"] for r in range(4)] for k in range(4)] for m in range(4)]
    A = [[ETA[m] * H[m][k] for k in range(4)] for m in range(4)]          # eta^{-1} H

    def mm(X, Y):
        return [[sum((X[i][k] * Y[k][j] for k in range(4)), R.zero) for j in range(4)] for i in range(4)]
    I4 = [[R.one if i == j else R.zero for j in range(4)] for i in range(4)]
    A2 = mm(A, A)
    # graded pieces: g^{-1} = (1 - eps A + eps^2 A^2) eta ;  entries (k: matrix) with [i][j] -> X[i][j]*ETA[j]
    ginv = {0: [[I4[i][j] * ETA[j] for j in range(4)] for i in range(4)],
            1: [[-A[i][j] * ETA[j] for j in range(4)] for i in range(4)],
            2: [[A2[i][j] * ETA[j] for j in range(4)] for i in range(4)]}
    t1 = sum((A[i][i] for i in range(4)), R.zero)
    t2 = sum((A2[i][i] for i in range(4)), R.zero)
    # sqrt(det(1+eps A)) = 1 + eps t1/2 + eps^2 (t1^2/8 - t2/4) + ...
    sq = {0: R.one, 1: t1 / 2, 2: t1 ** 2 / 8 - t2 / 4}
    # Gamma_k^a_{mn} (H-degree k): 1/2 ginv_k^{ab} (dH_{bn,m} + dH_{bm,n} - dH_{mn,b})
    G = {k: [[[sum((ginv[k][a][b] * (dH[b][n][m] + dH[b][m][n] - dH[m][n][b]) for b in range(4)), R.zero) / 2
               for n in range(4)] for m in range(4)] for a in range(4)] for k in range(3)}
    L = R.zero
    for i in range(3):
        for j in range(3 - i):
            for k in range(3 - i - j):
                l_ = 2 - i - j - k
                acc = R.zero
                for m in range(4):
                    for n in range(4):
                        gmn = ginv[j][m][n]
                        if not gmn:
                            continue
                        inner = R.zero
                        for a in range(4):
                            for b in range(4):
                                inner += G[k][a][m][b] * G[l_][b][n][a] - G[k][a][m][n] * G[l_][b][a][b]
                        acc += gmn * inner
                L += sq[i] * acc
    return sp.expand(L.as_expr())


def probe(e):
    phi, f = sp.Symbol("phi"), sp.symbols("f0:4")
    sub = {}
    for s in e.free_symbols:
        n = s.name
        if n.startswith("dH"):
            m, k, r = int(n[2]), int(n[3]), int(n.split("_")[1])
            sub[s] = ETA[m] * f[r] if m == k else 0
        elif n.startswith("H"):
            m, k = int(n[1]), int(n[2])
            sub[s] = ETA[m] * phi if m == k else 0
    return sp.expand(e.xreplace(sub)), phi, f


def main(export_path):
    t0 = time.time()
    exp_ = sp.sympify(Path(export_path).read_text())            # XVI-P p23 export (srepr)
    c4 = build()
    xo = None
    pr_c4, phi, f = probe(c4)
    xo = sp.expand(-3 * phi ** 2 * f[0] ** 2 + 3 * phi ** 2 * (f[1] ** 2 + f[2] ** 2 + f[3] ** 2))   # XVI-O probe
    s = sp.cancel(xo / pr_c4) if pr_c4 != 0 else None
    mine = sp.expand(s * c4)
    diff = sp.expand(mine - exp_)
    pr_exp, _, _ = probe(exp_)
    out = {"route": "independent Gamma-Gamma series (this file), not XVI-P quartic.py",
           "normalisation_from_XVI_O_probe": str(s), "n_terms_mine": len(sp.Add.make_args(mine)),
           "n_terms_export": len(sp.Add.make_args(exp_)), "full_polynomial_equal": diff == 0,
           "n_terms_difference": 0 if diff == 0 else len(sp.Add.make_args(diff)),
           "export_on_XVI_O_probe_equal": sp.expand(pr_exp - xo) == 0, "seconds": round(time.time() - t0, 1)}
    # mutation: one coefficient of the export changed must be detected
    t1 = sp.Add.make_args(exp_)[0]
    out["mutation_one_term_x2_detected"] = sp.expand(mine - (exp_ + t1)) != 0
    print(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    main(sys.argv[1])
