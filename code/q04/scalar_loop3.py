#!/usr/bin/env python3
"""XVI-Q Q04: scalar-loop kappa^3 graviton three-point UV pole (P.6 graphs G13 = V_h_phi_phi^3 triangle, G14 =
V_h_phi_phi + V_hh_phi_phi bubble; the V_hhh_phi_phi phi-phi tadpole is scaleless) and the matter-loop diffeomorphism
Ward identity at three points:
   Gamma3(delta0 h(q), e1(p1), e2(p2)) + Gamma2(delta1 h[xi, e1](q+p1), e2(p2)) + Gamma2(delta1 h[xi, e2](q+p2), e1(p1)) = 0,
   delta0 h = -i (q (x) xi + xi (x) q),  delta1 h[xi, e](k) = -i [ (k.xi) e + q (x) (e.xi) + (e.xi) (x) q ]   (L_xi g)
(W[g] of the scalar loop is diffeomorphism invariant; g = eta + h, kappa = 1).  Pole in units of U (uvpole.py).
Conventions tie: Gamma2 here equals I x XVI-P P18 D11 pole on all 14 structures (validated in q04 certs)."""
from __future__ import annotations
import sys
from pathlib import Path
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import uvpole as U  # noqa: E402

HALF = sp.Rational(1, 2)


def gamma2(eA, eB, k):
    """i Gamma2 pole / U: A carries incoming momentum k, B carries -k.  Scalar bubble, symmetry 1/2."""
    l = U.vec("l")
    num = U.V_hphiphi(eA, l, -(l + k)) * U.V_hphiphi(eB, l + k, -l)
    return U.loop_pole(num, [U.V(), k], 4, HALF)


def gamma3(legs):
    """legs = [(eps, k)] x 3 (incoming momenta, sum 0).  Triangles: both cyclic orders, factor 1/2;
    bubbles: X2 vertex with a pair + X1 vertex with the third, factor 1/2."""
    l = U.vec("l")
    tot = 0
    parts = {}
    for order in ((0, 1, 2), (0, 2, 1)):
        (eA, kA), (eB, kB), (eC, kC) = (legs[i] for i in order)
        num = U.V_hphiphi(eA, l, -(l + kA)) * U.V_hphiphi(eB, l + kA, -(l + kA + kB)) * U.V_hphiphi(eC, l + kA + kB, -l)
        v = U.loop_pole(num, [U.V(), kA, kA + kB], 4, HALF)
        parts[f"triangle{order}"] = v
        tot += v
    for c in range(3):
        a, b = [i for i in range(3) if i != c]
        (eC, kC), (ea, _), (eb, _) = legs[c], legs[a], legs[b]
        num = U.V_hphiphi(eC, l, -(l + kC)) * U.V_hhphiphi(ea, eb, l + kC, -l)
        v = U.loop_pole(num, [U.V(), kC], 4, HALF)
        parts[f"bubble_X2{(a, b)}_X1{c}"] = v
        tot += v
    return sp.expand(tot), parts


def pv(*pairs):
    """polarisation-side vector WITHOUT the Taylor parameter (constant w.r.t. the loop expansion)."""
    return U.V(dict(pairs))


def ward(drop_bubbles=False, lie_sign=1):
    p1, p2 = U.vec("p1"), U.vec("p2")
    q = -(p1 + p2)
    qc = pv(("p1", -1), ("p2", -1))
    xi = pv(("xi", 1))
    e1, e2 = U.rank1("a1"), U.rank1("a2")
    d0 = U.Eps([(-2 * sp.I, qc, xi)])
    g3, parts = gamma3([(d0, q), (e1, p1), (e2, p2)])
    if drop_bubbles:
        g3 = sp.expand(sum(v for k, v in parts.items() if k.startswith("triangle")))

    def d1(a, kvec_c):
        ax = U.dot(U.vec(a), xi)                          # (e.xi) = (a.xi) a  (4S)
        return U.Eps([(-sp.I * U.dot(kvec_c, xi), U.vec(a), U.vec(a)), (-2 * sp.I * ax * lie_sign, qc, U.vec(a))])
    w = g3 + gamma2(d1("a1", pv(("p1", 1))), e2, -p2) + gamma2(d1("a2", pv(("p2", 1))), e1, -p1)
    return sp.expand(w), g3


if __name__ == "__main__":
    import time
    t0 = time.time()
    w, g3 = ward()
    print("ward residual (pole/U):", sp.factor(w) if w != 0 else 0, time.time() - t0)
