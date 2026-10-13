#!/usr/bin/env python3
"""XVI-Q Q04 independent one-loop UV-pole engine (own Feynman rules; NOT XVI-P exports).

Prescription (XVI-P declared, DRED): loop momentum l in QdS; external momenta p_i in QdS; vertex (graviton/ghost)
indices in Q4S = QdS (+) Q(4-d)S.  Polarisation vectors a carry a hat part (QdS) and a tilde part (evanescent).
Invariants: H(x,y) = hat dot, T(x,y) = tilde dot (only between polarisation vectors); 4S dot = H + T; momenta are hat.
d is the XVI-P opaque dimension (qguard); it appears only through symmetric integration
    < l^{mu1} ... l^{mu2n} > = sym(ghat)^n (l^2)^n / (d (d+2) ... (d+2n-2)).
UV pole: Taylor expansion in the external momenta (scaling parameter t) to the superficial degree omega (exact
algebraic identity per term); the t^omega coefficient is homogeneous of l-degree -4 and its UV pole is
<...> * U,  U := pole of int d^dl/(2pi)^d (l^2+i0)^-2 = i/(16 pi^2 eps).  Terms of other l-degree are power divergent
(scaleless after IR rearrangement: pole proportional to a regulator mass, vanishing at M = 0) or UV finite.
Massless: no mass anywhere (XVI-O O.1 m = 0)."""
from __future__ import annotations
import itertools, sys  # noqa: E401
from functools import lru_cache
from pathlib import Path
import sympy as sp

Q = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Q / "code" / "q01"))
import qguard  # noqa: E402,F401   (installs the d guard)
d = qguard.d
t = sp.Symbol("t")
L = sp.Symbol("L")           # l.l
MOM = ("p1", "p2")           # independent external momenta (QdS)


def is_mom(n):
    return n == "l" or n in MOM


def H(x, y):
    x, y = sorted((x, y))
    return sp.Symbol(f"H_{x}_{y}")


def T(x, y):
    x, y = sorted((x, y))
    return sp.Symbol(f"T_{x}_{y}")


def S(x):                    # l.x
    return sp.Symbol(f"s_{x}")


class V(dict):
    """vector = {basis name: coefficient}; basis 'l', 'p1', 'p2' (hat) or polarisation names (hat+tilde)."""
    def __add__(self, o):
        r = V(self)
        for k, v in o.items():
            r[k] = r.get(k, 0) + v
        return r

    def __neg__(self):
        return V({k: -v for k, v in self.items()})

    def __sub__(self, o):
        return self + (-o)

    def __rmul__(self, c):
        return V({k: c * v for k, v in self.items()})


def vec(n, c=1):
    return V({n: (t * c if n in MOM else c)})      # momenta carry the Taylor parameter t


def bdot(x, y, hat_only=False):
    if x == "l" and y == "l":
        return L
    if x == "l" or y == "l":
        o = y if x == "l" else x
        return S(o)
    if is_mom(x) or is_mom(y) or hat_only:
        return H(x, y)
    return H(x, y) + T(x, y)


def dot(u, w, hat_only=False):
    return sp.expand(sum(cu * cw * bdot(a, b, hat_only) for a, cu in u.items() for b, cw in w.items()))


class Eps(list):
    """symmetric 4S tensor = sum c * (u (x) w + w (x) u)/2 ; items (c, u, w)."""
    def __call__(self, x, y):                       # eps(x, y)
        return sp.expand(sum(sp.Rational(1, 2) * c * (dot(u, x) * dot(w, y) + dot(w, x) * dot(u, y)) for c, u, w in self))

    def tr(self):
        return sp.expand(sum(c * dot(u, w) for c, u, w in self))

    def vdot(self, x):                              # (eps . x)_mu as a list of (coef, vector)
        return [(sp.Rational(1, 2) * c * dot(w, x), u) for c, u, w in self] + [(sp.Rational(1, 2) * c * dot(u, x), w) for c, u, w in self]


def rank1(a):
    return Eps([(1, vec(a), vec(a))])


def pair_contract(e1, e2, x, y):
    """(e1 . e2)(x, y) = e1_{x rho} e2^{rho}_y (4S contraction over rho)."""
    tot = 0
    for c1, v1 in e1.vdot(x):
        for c2, v2 in e2.vdot(y):
            tot += c1 * c2 * dot(v1, v2)
    return sp.expand(tot)


def tr2(e1, e2):
    """tr4(e1 e2)."""
    tot = 0
    for c1, u1, w1 in e1:
        for c2, u2, w2 in e2:
            tot += sp.Rational(1, 2) * c1 * c2 * (dot(u1, u2) * dot(w1, w2) + dot(u1, w2) * dot(w1, u2))
    return sp.expand(tot)


# ---------------- scalar sector:  L_phi = (1/2) sqrt(-g) g^{mn} d_m phi d_n phi,  g = eta + kappa h  (eta + - - -)
#   X^{mn} = sqrt(-g) g^{mn} = eta + kappa X1 + kappa^2 X2 + ... ; X1 = (1/2) tr h eta - h ; X2 = h.h - (1/2) tr h h + (tr h^2/8 - tr(h.h)/4) eta
TRACE_COEF = sp.Rational(1, 2)    # mutation hook (q04 certs): sqrt(-g) -> 1 + TRACE_COEF tr h


def X1(e, x, y):
    return sp.expand(TRACE_COEF * e.tr() * dot(x, y) - e(x, y))


def X2cross(e1, e2, x, y):
    """coefficient of t1 t2 in X2[t1 e1 + t2 e2](x, y)."""
    return sp.expand(pair_contract(e1, e2, x, y) + pair_contract(e2, e1, x, y)
                     - sp.Rational(1, 2) * (e1.tr() * e2(x, y) + e2.tr() * e1(x, y))
                     + (sp.Rational(1, 4) * e1.tr() * e2.tr() - sp.Rational(1, 2) * tr2(e1, e2)) * dot(x, y))


def X3(es, x, y):
    """fully polarised cubic term of sqrt(-g) g^{mn} (coefficient of t1 t2 t3) -- used for the V_hhh_phi_phi vertex."""
    raise NotImplementedError


def V_hphiphi(e, k1, k2):
    """vertex i * delta^3 S / (delta h delta phi delta phi), scalar momenta k1, k2 incoming: -i X1(k1, k2) (kappa^1)."""
    return -sp.I * X1(e, k1, k2)


def V_hhphiphi(e1, e2, k1, k2):
    return -sp.I * X2cross(e1, e2, k1, k2)


# ---------------- symmetric integration and UV pole
@lru_cache(maxsize=None)
def _avg_coeff(key):
    """<prod_j (l.v_j)^{k_j}> / (l^2)^{K/2} as a polynomial in hat dots; key = tuple((name, k_j))."""
    names = [n for n, _ in key]
    ks = [k for _, k in key]
    K = sum(ks)
    if K % 2:
        return sp.Integer(0)
    tt = sp.symbols(f"tt0:{len(names)}")
    x2 = sp.expand(sum(tt[i] * tt[j] * H(names[i], names[j]) for i in range(len(names)) for j in range(len(names))))
    P = sp.Poly(sp.expand(x2 ** (K // 2)), *tt)
    c = P.coeff_monomial(tuple(ks))
    num = sp.Integer(1)
    for k in ks:
        num *= sp.factorial(k)
    dfac = sp.Integer(1)
    for i in range(K // 2):
        dfac *= (d + 2 * i)
    return sp.expand(c * num / sp.factorial(K) * sp.factorial2(K - 1)) / dfac


def integrate_pole(expr):
    """expr: polynomial in s_x and L^{+-1} (homogeneous l-degree -4 terms selected); returns pole / U."""
    expr = sp.expand(expr)
    svars = sorted({s for s in expr.free_symbols if s.name.startswith("s_")}, key=lambda s: s.name)
    out = 0
    for term in sp.Add.make_args(expr):
        powd = term.as_powers_dict()
        kL = powd.get(L, 0)
        ks = [(s.name[2:], int(powd.get(s, 0))) for s in svars if powd.get(s, 0)]
        K = sum(k for _, k in ks)
        if K + 2 * kL != -4:
            continue
        coef = term / (L ** kL)
        for s in svars:
            coef = coef / s ** powd.get(s, 0)
        out += coef * (_avg_coeff(tuple(ks)) if ks else 1)
    return out


def taylor_denominators(dens, omega):
    """prod_j 1/(l+K_j)^2 expanded in t to order omega: sum_n (-(2 l.K + K^2)/L)^n / L."""
    tot = sp.Integer(1)
    for K in dens:
        lk = dot(vec("l"), K)
        kk = dot(K, K)
        x = sp.expand(2 * lk + kk)
        ser = sum((-x) ** n / L ** (n + 1) for n in range(omega + 1))
        tot = sp.expand(tot * ser)
        tot = trunc(tot, omega)
    return tot


def trunc(e, n):
    e = sp.expand(e)
    return sum(term for term in sp.Add.make_args(e) if sp.degree(term, t) <= n) if e != 0 else e


def loop_pole(numer, dens, omega, prefactor):
    """UV pole / U of prefactor * int numer * prod i/(l+K)^2 ; numer contains the vertex factors."""
    num = trunc(sp.expand(numer), omega)
    den = taylor_denominators(dens, omega)
    full = sp.expand(num * den)
    top = sp.expand(full).coeff(t, omega)
    return sp.expand(prefactor * sp.I ** len(dens) * integrate_pole(top))
