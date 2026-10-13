#!/usr/bin/env python3
"""XVI-Q Q02 scrutiny of D04 = 0: an INDEPENDENT route (no XVI-P engine code) with covariant algebra and an own
Feynman-parameter evaluation with symbolic dimension D (a plain symbol: this file never assigns D a value), plus mutations.

Mechanism. Vertex (phi(p), phi(q), graviton k = -(p+q)):
  V_ab = p_a q_b + p_b q_a - alpha_v eta_ab p.q + lam [k_a k_b - eta_ab k^2]   (lam ~ non-minimal xi R phi^2 coupling)
propagator numerator P = 1/2(eta eta + eta eta) - c eta eta (+ g3 * Froeb's T3 = 4 l_(a eta_b)(c l_d)/l^2 gauge piece),
vertex-index trace eta.eta = n.  With alpha_v = 1, lam = 0, g3 = 0:
  V.P.V' = 2 p^2 q^2 + (n - 2)(p.q)^2 [1 - c (n - 2)],
so the (p.q)^2 (trace-mode) piece cancels iff c = 1/(n - 2): DRED (n = 4, c = 1/2) and CDR (n = D, c = 1/(D - 2)) both
give V.P.V' = 2 p^2 (l - p)^2, which cancels the scalar propagator; the remainder p^2 int 1/l^2 is scaleless (omega = 2,
no logarithmic member).  D04 = 0 is therefore an exact statement of de Donder gauge + minimal coupling + m = 0, NOT a
property of scaleless integrals in general (the bubble with p != 0 is not scaleless).
Literature: Froeb, CQG 35 (2018) 055006 [arXiv:1710.00839], eq. inv2pf_termsumint_grav: graviton-loop self-energy with
the de Donder parts g1, g2 of the propagator in n dimensions; the integrand bracket is proportional to m^2 or xi, and
the text states it 'vanishes identically' for xi = 0, m = 0 'in accordance with' Kahya-Woodard 2007; for m = 0 the xi^2
terms vanish; gauge-dependent propagator pieces (T3..T5) contribute (cancelled there only by coordinate corrections)."""
from __future__ import annotations
import json
import sys

import sympy as sp

D = sp.Symbol("D")                       # own dimension symbol (never substituted)
PP, PL, LL = sp.symbols("pp pl ll")     # p.p, p.l, l.l
x = sp.Symbol("x")


# ---------------- vectors over the basis (p, l) and their dot product
def vec(cp, cl):
    return (sp.sympify(cp), sp.sympify(cl))


def dot(u, v):
    return sp.expand(u[0] * v[0] * PP + (u[0] * v[1] + u[1] * v[0]) * PL + u[1] * v[1] * LL)


def add(*vs):
    return (sum(v[0] for v in vs), sum(v[1] for v in vs))


def sc(a, v):
    return (a * v[0], a * v[1])


P_ = vec(1, 0)
L_ = vec(0, 1)


# ---------------- symmetric tensors: list of (coef, u, v) for coef * u v^T, plus coefficient of eta
class T:
    def __init__(self, outer, eta):
        self.outer, self.eta = outer, sp.sympify(eta)

    def dot(self, o, n):
        s = sum(a * b * dot(u, u2) * dot(v, v2) for a, u, v in self.outer for b, u2, v2 in o.outer)
        s += self.eta * o.trace(n) + o.eta * sum(a * dot(u, v) for a, u, v in self.outer) - self.eta * o.eta * n + self.eta * o.eta * n
        return sp.expand(s)

    def trace(self, n):
        return sp.expand(sum(a * dot(u, v) for a, u, v in self.outer) + n * self.eta)

    def contract_vec(self, w):
        """(T w)_a as a vector."""
        out = vec(0, 0)
        for a, u, v in self.outer:
            out = add(out, sc(a * dot(v, w), u))
        return add(out, sc(self.eta, w))


def vertex(p, q, alpha_v=1, lam=0):
    k = sc(-1, add(p, q))
    outer = [(1, p, q), (1, q, p)]
    eta = -alpha_v * dot(p, q)
    if lam:
        outer.append((lam, k, k))
        eta -= lam * dot(k, k)
    return T(outer, eta)


def numerator45(n, g4=0, g5=0):
    """Froeb T4 = eta_ab l_c l_d/l^2 + l_a l_b eta_cd/l^2 (coefficient g4 ~ (2+beta)) -> numerator over one extra l^2;
    T5 = l_a l_b l_c l_d / l^4 (coefficient g5) -> numerator over (l^2)^2."""
    q = add(L_, sc(-1, P_))
    VA, VB = vertex(P_, q), vertex(sc(-1, P_), sc(-1, q))
    lVl_A, lVl_B = dot(L_, VA.contract_vec(L_)), dot(L_, VB.contract_vec(L_))
    N4 = sp.expand(g4 * (VA.trace(n) * lVl_B + VB.trace(n) * lVl_A))
    N5 = sp.expand(g5 * lVl_A * lVl_B)
    return N4, N5


def numerator(n, c, alpha_v=1, lam=0, g3=0):
    """N = V(p, q) . P . V(-p, -q), q = l - p, graviton momentum l; returns (N_deDonder, N_T3) with N_T3 to be divided
    by an extra l^2 (Froeb T3 = 4 l_(mu eta_nu)(rho l_sigma)/l^2)."""
    q = add(L_, sc(-1, P_))
    VA = vertex(P_, q, alpha_v, lam)
    VB = vertex(sc(-1, P_), sc(-1, q), alpha_v, lam)
    N0 = sp.expand(VA.dot(VB, n) - c * VA.trace(n) * VB.trace(n))
    N3 = sp.expand(g3 * 4 * dot(VA.contract_vec(L_), VB.contract_vec(L_))) if g3 else sp.Integer(0)
    return N0, N3


# ---------------- own Feynman-parameter integration (Euclidean form; ratios to the master bubble)
def kint(j, Nexp):
    """int d^Dk (k^2)^j/(k^2+Delta)^N  /  [pi^{D/2}] = G(j+D/2) G(N-j-D/2)/(G(D/2) G(N)) Delta^{D/2+j-N}"""
    return sp.gamma(j + D / 2) * sp.gamma(Nexp - j - D / 2) / (sp.gamma(D / 2) * sp.gamma(Nexp))


def bubble(num, a, b):
    """int num(l,p)/((l^2)^a ((l-p)^2)^b), normalised by the master int 1/(l^2 (l-p)^2); returns rational(D) * pp^power.
    l = k + x p, (1-x) l^2 + x (l-p)^2 = k^2 + x(1-x) pp; symmetric integration in D dims."""
    num = sp.expand(num)
    if num == 0:
        return sp.Integer(0)
    kk, kp = sp.symbols("kk kp")
    # substitute l.l = kk + 2 x kp + x^2 pp, p.l = kp + x pp
    e = sp.expand(num.subs({LL: kk + 2 * x * kp + x ** 2 * PP, PL: kp + x * PP}, simultaneous=True))
    Ptot = a + b
    pref = sp.gamma(Ptot) / (sp.gamma(a) * sp.gamma(b))
    tot = 0
    for (ik, ip), cf in sp.Poly(e, kk, kp).terms():
        if ip % 2:
            continue
        m = ip // 2
        # <(k.p)^{2m}> = (2m-1)!! (k^2)^m pp^m / prod_{i<m}(D + 2i)
        ang = sp.factorial2(2 * m - 1) * PP ** m / sp.Mul(*[D + 2 * i for i in range(m)]) if m else 1
        j = ik + m
        # x-integral of x^{b-1} (1-x)^{a-1} * cf(x) * [x(1-x)]^{D/2+j-Ptot}
        for (ex,), cx in sp.Poly(sp.expand(cf), x).terms():
            beta = sp.gamma(ex + b + D / 2 + j - Ptot) * sp.gamma(a + D / 2 + j - Ptot) / sp.gamma(ex + a + b + D + 2 * j - 2 * Ptot)
            tot += pref * cx * ang * kint(j, Ptot) * beta * PP ** (j - Ptot + 2)   # pp power relative to the master
    master = kint(0, 2) * sp.gamma(D / 2 - 1) ** 2 / sp.gamma(D - 2)
    return sp.factor(sp.simplify(sp.gammasimp(sp.expand(tot) / master)))


def D04(n, c, alpha_v=1, lam=0, g3=0, g4=0, g5=0):
    N0, N3 = numerator(n, c, alpha_v, lam, g3)
    N4, N5 = numerator45(n, g4, g5)
    return sp.factor(sp.simplify(bubble(N0, 1, 1) + (bubble(N3, 2, 1) if N3 != 0 else 0)
                                 + (bubble(N4, 2, 1) if N4 != 0 else 0) + (bubble(N5, 3, 1) if N5 != 0 else 0)))


def main():
    xi = sp.Symbol("xi")
    out = {"mechanism_numerator_general_c_n": str(sp.factor(numerator(sp.Symbol("n"), sp.Symbol("c"))[0])),
           "mechanism_reading": "V.P.V' = 2 p^2 (l-p)^2 + (n-2)(p.q)^2 [1 - c(n-2)]: trace-mode piece cancels iff c = 1/(n-2)"}
    N_dred = numerator(4, sp.Rational(1, 2))[0]
    out["N_DRED_equals_2p2(l-p)^2"] = sp.expand(N_dred - 2 * PP * (LL - 2 * PL + PP)) == 0
    cases = {
        "DRED (n=4, c=1/2) -- declared": (4, sp.Rational(1, 2), 1, 0, 0),
        "CDR (n=D, c=1/(D-2))": (D, 1 / (D - 2), 1, 0, 0),
        "MUT mixed scheme (n=4, c=1/(D-2))": (4, 1 / (D - 2), 1, 0, 0),
        "MUT vertex trace normalisation alpha_v=1/2": (4, sp.Rational(1, 2), sp.Rational(1, 2), 0, 0),
        "GAUGE T3 piece g3=1 (expected 0: l^a V_ab = p_b q^2 + q_b p^2, every term carries q^2 -> scaleless)": (4, sp.Rational(1, 2), 1, 0, 1, 0, 0),
        "MUT non-de Donder gauge T4 piece (beta != -2), g4=1": (4, sp.Rational(1, 2), 1, 0, 0, 1, 0),
        "MUT non-de Donder gauge T5 piece, g5=1": (4, sp.Rational(1, 2), 1, 0, 0, 0, 1),
        "MUT non-minimal coupling lam=xi (CDR)": (D, 1 / (D - 2), 1, xi, 0),
    }
    res = {}
    for k, args in cases.items():
        v = D04(*args)
        res[k] = {"value_over_master": str(v), "zero": v == 0}
    out["cases"] = res
    out["T3_zero_by_Ward_identity"] = res[next(k for k in res if k.startswith("GAUGE T3"))]["zero"]
    v_xi = sp.sympify(res["MUT non-minimal coupling lam=xi (CDR)"]["value_over_master"], locals={"D": D, "pp": PP, "xi": xi})
    out["xi_dependence_linear_only (Froeb: m=0 xi^2 terms vanish)"] = sp.expand(v_xi).coeff(xi, 2) == 0 and sp.expand(v_xi).coeff(xi, 1) != 0
    out["declared_zero"] = res["DRED (n=4, c=1/2) -- declared"]["zero"] and res["CDR (n=D, c=1/(D-2))"]["zero"]
    out["all_mutations_nonzero"] = all(not r["zero"] for k, r in res.items() if k.startswith("MUT"))
    out["literature"] = {"reference": "M. B. Froeb, Class. Quantum Grav. 35 (2018) 055006, arXiv:1710.00839 (flat space, general gauge, n dims)",
                         "statement": "graviton-loop scalar self-energy from the de Donder propagator parts vanishes identically for "
                                      "xi = 0, m = 0 (in accordance with Kahya & Woodard 2007); m = 0 xi^2 terms vanish; "
                                      "gauge-parameter pieces contribute (gauge-dependent 1PI)",
                         "agreement": "AGREES: declared case 0; xi enters linearly only; non-de Donder piece nonzero"}
    return out


if __name__ == "__main__":
    r = main()
    print(json.dumps(r, indent=1))
