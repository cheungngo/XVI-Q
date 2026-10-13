#!/usr/bin/env python3
"""XVI-Q Q06: quartic Einstein-Hilbert vertex V4_grav from the XVI-P action (same conventions as p6q05/gravity.py):
S_EH = -(2/kappa^2) int L_GG, L_GG = sqrt(-g) g^{mn} g^{ac} g^{bd} (G_{c,mb} G_{d,na} - G_{c,mn} G_{d,ab}), g = eta + kappa h.
The prefactor P = sqrt(-g) g^{mn} g^{ac} g^{bd} is expanded to second order (polarised coefficient of t1 t2):
  sqrt(-g):  -1/2 tr(E1 E2) + 1/4 trE1 trE2 ;  g^{-1}: -E (1st), E1 E2 + E2 E1 (2nd).
Validation (no non-scaleless XVI-P graph contains V4): the Einstein-Hilbert gauge identity order by order,
  delta0 S_{n+1} + delta1 S_n = 0,  delta0 h = d xi + d xi,  delta1 h = L_xi h,
checked exactly for n = 2 (V3 vs V2, machinery control) and n = 3 (V4 vs V3), plus total permutation symmetry of V4."""
from __future__ import annotations
import itertools, sys  # noqa: E401
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "p6q05"))
import gravity as G  # noqa: E402
import uvpole as U  # noqa: E402
import tensor as TT  # noqa: E402
from tensor import ETA, ev  # noqa: E402

I = sp.I
L = G.L


def _slot(M, x, y):
    return [("m", M, x, y)]


def _chain(M1, M2, x, y):
    z = L()
    return [("m", M1, x, z), ("m", M2, z, y)]


def prefactor_terms2(E1, E2, m, n, a, c, b, d):
    """coefficient of t1 t2 in sqrt(-g) g^{mn} g^{ac} g^{bd} at h = t1 E1 + t2 E2."""
    slots = [(m, n), (a, c), (b, d)]
    i1, i2 = L(), L()
    s2 = sp.Rational(1, 4) * E1.tr() * E2.tr() - sp.Rational(1, 2) * ev((1, [("m", E1, i1, i2), ("m", E2, i2, i1)]))
    out = []

    def build(fill):
        f = []
        for k, (x, y) in enumerate(slots):
            f += fill.get(k, None) or _slot(ETA(), x, y)
        return f
    out.append((s2, build({})))
    for k, (x, y) in enumerate(slots):
        out.append((-sp.Rational(1, 2) * E1.tr(), build({k: _slot(E2, x, y)})))
        out.append((-sp.Rational(1, 2) * E2.tr(), build({k: _slot(E1, x, y)})))
        out.append((1, build({k: _chain(E1, E2, x, y)})))
        out.append((1, build({k: _chain(E2, E1, x, y)})))
    for k1, k2 in itertools.combinations(range(3), 2):
        (x1, y1), (x2, y2) = slots[k1], slots[k2]
        out.append((1, build({k1: _slot(E1, x1, y1), k2: _slot(E2, x2, y2)})))
        out.append((1, build({k1: _slot(E2, x1, y1), k2: _slot(E1, x2, y2)})))
    return out


def LGG_piece2(pref, w1, w2):
    m, n, a, b, c, d = (L() for _ in range(6))
    P = prefactor_terms2(pref[0], pref[1], m, n, a, c, b, d)
    t1 = G._mul(G._mul(P, G.gamma_terms(w1[0], w1[1], c, m, b)), G.gamma_terms(w2[0], w2[1], d, n, a))
    m2, n2, a2, b2, c2, d2 = (L() for _ in range(6))
    P2 = prefactor_terms2(pref[0], pref[1], m2, n2, a2, c2, b2, d2)
    t2 = G._mul(G._mul(P2, G.gamma_terms(w1[0], w1[1], c2, m2, n2)), G.gamma_terms(w2[0], w2[1], d2, a2, b2))
    return ev(*t1) - ev(*t2)


def LGG4(w1, w2, w3, w4):
    ws = [w1, w2, w3, w4]
    tot = 0
    for i, j in itertools.permutations(range(4), 2):
        rest = [ws[k][0] for k in range(4) if k not in (i, j)]
        tot += LGG_piece2(rest, ws[i], ws[j])
    return sp.expand(tot)


def V4_grav(w1, w2, w3, w4):
    """i x delta^4 S_EH (kappa^2 stripped): S_EH quartic = -2 kappa^2 L_GG^(4)."""
    return sp.expand(I * (-2) * LGG4(w1, w2, w3, w4))


def V2_EH(w1, w2):
    return sp.expand(I * (-2) * G.LGG2(w1, w2))


def gauge_check(n):
    """delta0 V_{n+1} + delta1 V_n with xi (momentum q) and n graviton waves a_i (momenta k_i, sum k + q = 0)."""
    xi, q = U.vec("xi"), U.vec("q")
    ks = [U.vec("k1"), U.vec("k2")] + ([U.vec("k3")] if n == 3 else [])
    ks = ks[:n - 1] + [-(sum(ks[:n - 1], U.V()) + q)]
    Es = [TT.sym(U.vec(f"a{i+1}"), U.vec(f"a{i+1}")) for i in range(n)]
    W = [(Es[i], ks[i]) for i in range(n)]
    lin = (TT.sym(q, xi, -2 * I), q)
    Vn1 = G.V3_grav if n == 2 else V4_grav
    Vn = V2_EH if n == 2 else G.V3_grav
    tot = Vn1(lin, *W)
    for i in range(n):
        Wi = list(W)
        Wi[i] = (G.Lc_h(Es[i], ks[i], xi, q), ks[i] + q)
        tot += Vn(*Wi)
    return sp.expand(tot)


def symmetry_check():
    ks = [U.vec("k1"), U.vec("k2"), U.vec("k3")]
    ks.append(-(ks[0] + ks[1] + ks[2]))
    W = [(TT.sym(U.vec(f"a{i+1}"), U.vec(f"a{i+1}")), ks[i]) for i in range(4)]
    base = V4_grav(*W)
    sw = V4_grav(W[1], W[0], W[3], W[2])
    return sp.expand(base - sw) == 0, len(sp.Add.make_args(base))


if __name__ == "__main__":
    import json, time  # noqa: E401
    t0 = time.time()
    r = {"gauge_n2_V3_vs_V2_zero": gauge_check(2) == 0}
    print(r, time.time() - t0, flush=True)
    g3 = gauge_check(3)
    r["gauge_n3_V4_vs_V3_zero"] = g3 == 0
    r["gauge_n3_residual_terms"] = len(sp.Add.make_args(g3))
    print(r, time.time() - t0, flush=True)
    r["V4_permutation_symmetric"], r["V4_terms_generic"] = symmetry_check()
    r["seconds"] = round(time.time() - t0, 1)
    print(r, flush=True)
    (HERE.parents[1] / "build" / "logs" / "q06_V4_validation.json").write_text(json.dumps(r, indent=1) + "\n")
