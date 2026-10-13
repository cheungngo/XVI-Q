#!/usr/bin/env python3
"""XVI-Q Q05 own gravity/ghost Feynman rules from the XVI-P declared action (conventions.py, read-only):
  S_EH = -(2/kappa^2) int L_GG,  L_GG = sqrt(-g) g^{mn}(G^a_{mb} G^b_{na} - G^a_{mn} G^b_{ab}),  g = eta + kappa h (+---)
  S_gf (b eliminated, XVI-P P13 'b(2D-b): eliminating b gives +D_mu D^mu'):  int D_mu[h] D^mu[h],
        D^mu[h] = d_nu h^{mu nu} - (1/2) d^mu tr h
  S_gh  = int i cbar_mu D^mu[ (d c + d c) + kappa L_c h ]   (XVI-P: kinetic 'i Box', vertex '2 i a D[L_c h]', a = 1/2)
  S_K   = int K^{mn} (1/2) [ (d_m c_n + d_n c_m) + kappa (L_c h)_{mn} ]   (K s h, s h = (1/2) L_c(eta + kappa h))
  scalar: (1/2) sqrt(-g) g^{mn} d phi d phi.
Plane waves: field_i = pol_i exp(-i k_i x) (all momenta incoming), d -> -i k.  A vertex = i x (multilinear
coefficient of the action).  kappa powers are tracked by the caller (vertex functions are at kappa = 1).
Internal lines are closed by PROBES: the two ends of a line get probe polarisations; the product of vertices is
contracted with the propagator numerator (graviton: i/(Z k^2) [ (eta eta + eta eta)/2 - (1/2) eta eta ], Q4S indices;
ghost: <c cbar> = i K^{-1}; scalar: i/k^2).  Normalisations Z etc. are DERIVED in derive_propagators()."""
from __future__ import annotations
import itertools, re, sys  # noqa: E401
from pathlib import Path
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "q04"))
import uvpole as U  # noqa: E402
import tensor as TT  # noqa: E402
from tensor import Mat, ETA, ev  # noqa: E402

I = sp.I
_lab = itertools.count()


def L():
    return f"i{next(_lab)}"


# ------------------------------------------------------------------ Christoffel (first kind), one plane wave
def gamma_terms(E, k, c, m, b):
    """Gamma_{c,mb} = 1/2 (d_m h_cb + d_b h_cm - d_c h_mb), d -> -i k: list of (coef, factors)."""
    h = -I / 2
    return [(h, [("v", k, m), ("m", E, c, b)]), (h, [("v", k, b), ("m", E, c, m)]), (-h, [("v", k, c), ("m", E, m, b)])]


def _mul(t1, t2):
    return [(a * b, fa + fb) for a, fa in t1 for b, fb in t2]


def prefactor_terms(order_wave, m, n, a, c, b, d):
    """sqrt(-g) g^{mn} g^{ac} g^{bd}: order 0 (order_wave None) or first order in wave (E, k)."""
    if order_wave is None:
        return [(1, [("m", ETA(), m, n), ("m", ETA(), a, c), ("m", ETA(), b, d)])]
    E = order_wave
    trE = E.tr()
    return [(trE / 2, [("m", ETA(), m, n), ("m", ETA(), a, c), ("m", ETA(), b, d)]),
            (-1, [("m", E, m, n), ("m", ETA(), a, c), ("m", ETA(), b, d)]),
            (-1, [("m", ETA(), m, n), ("m", E, a, c), ("m", ETA(), b, d)]),
            (-1, [("m", ETA(), m, n), ("m", ETA(), a, c), ("m", E, b, d)])]


def LGG_piece(pref_wave, w1, w2):
    """L_GG with the prefactor at order 0/1 (pref_wave) and Gamma's from waves w1, w2 (each (E, k))."""
    m, n, a, b, c, d = (L() for _ in range(6))
    P = prefactor_terms(pref_wave, m, n, a, c, b, d)
    t1 = _mul(_mul(P, gamma_terms(w1[0], w1[1], c, m, b)), gamma_terms(w2[0], w2[1], d, n, a))
    m2, n2, a2, b2, c2, d2 = (L() for _ in range(6))
    P2 = prefactor_terms(pref_wave, m2, n2, a2, c2, b2, d2)
    t2 = _mul(_mul(P2, gamma_terms(w1[0], w1[1], c2, m2, n2)), gamma_terms(w2[0], w2[1], d2, a2, b2))
    return ev(*t1) - ev(*t2)


def LGG2(w1, w2):
    """multilinear coefficient of t1 t2 in L_GG (both orders of the two Gammas)."""
    return sp.expand(LGG_piece(None, w1, w2) + LGG_piece(None, w2, w1))


def LGG3(w1, w2, w3):
    ws = [w1, w2, w3]
    tot = 0
    for p in itertools.permutations(range(3)):
        tot += LGG_piece(ws[p[0]][0], ws[p[1]], ws[p[2]])
    return sp.expand(tot)


def V3_grav(w1, w2, w3):
    """i x delta^3 S_EH (kappa^1 stripped): S_EH cubic = -2 kappa L_GG^(3)."""
    return sp.expand(I * (-2) * LGG3(w1, w2, w3))


def Dvec(E, k):
    """D^mu[h] for one wave as a vector: -i (E.k - (1/2) k tr E)."""
    return (-I) * (E.left(k) + (-sp.Rational(1, 2) * E.tr()) * k)


def S2_EH_gf(E1, E2, k):
    """delta^2 (S_EH + S_gf) with wave 1 (E1, k), wave 2 (E2, -k), kappa^0."""
    eh = -2 * LGG2((E1, k), (E2, -k))
    gf = 2 * U.dot(Dvec(E1, k), Dvec(E2, -k))
    return sp.expand(eh), sp.expand(gf)


# ------------------------------------------------------------------ scalar vertices on Mat polarisations
def X1(E, x, y):
    return sp.expand(sp.Rational(1, 2) * E.tr() * U.dot(x, y) - U.dot(E.left(x), y))


def X2cross(E1, E2, x, y):
    i1, i2 = L(), L()
    pc = lambda A, B: ev((1, [("v", x, i1), ("m", A, i1, i2), ("m", B, i2, L_ := L()), ("v", y, L_)]))  # noqa: E731
    return sp.expand(pc(E1, E2) + pc(E2, E1) - sp.Rational(1, 2) * (E1.tr() * U.dot(E2.left(x), y) + E2.tr() * U.dot(E1.left(x), y))
                     + (sp.Rational(1, 4) * E1.tr() * E2.tr() - sp.Rational(1, 2) * ev((1, [("m", E1, i1, i2), ("m", E2, i2, i1)]))) * U.dot(x, y))


def V_hphiphi(E, k1, k2):
    return -I * X1(E, k1, k2)


def V_hhphiphi(E1, E2, k1, k2):
    return -I * X2cross(E1, E2, k1, k2)


# ------------------------------------------------------------------ ghost / antifield vertices
def Lc_h(E, k, cvec, q):
    """(L_c h)_{mn} as a Mat for wave (E, k) and ghost polarisation cvec with momentum q:
    c^r d_r h_mn + h_rn d_m c^r + h_mr d_n c^r,  d -> -i (k on h, q on c)."""
    out = Mat(0, [(c0 * (-I) * U.dot(cvec, k), u, w) for c0, u, w in E.items])
    if E.eta != 0:
        out.eta = E.eta * (-I) * U.dot(cvec, k)
    Ec = E.left(cvec)                                      # (c.E)_n  (E symmetric)
    out.items += [(-I, q, Ec), (-I, Ec, q)]
    return out


def V_cb_c_h(cbv, cv, q, E, k):
    """i x coefficient of cbar(cbv) c(cv, q) h(E, k) in S_gh (kappa^1): i * i cbar.D[L_c h] with P = q + k."""
    M = Lc_h(E, k, cv, q)
    P = q + k
    D = (-I) * (M.left(P) + (-sp.Rational(1, 2) * M.tr()) * P)
    return sp.expand(I * I * U.dot(cbv, D))


def V_K_c_h(KE, cv, q, E, k):
    """i x coefficient of K(KE) c(cv, q) h(E, k) in S_K: i (1/2) K^{mn} (L_c h)_{mn} (K a symmetric 'polarisation' Mat)."""
    M = Lc_h(E, k, cv, q)
    i1, i2 = L(), L()
    return sp.expand(I * sp.Rational(1, 2) * ev((1, [("m", KE, i1, i2), ("m", M, i2, i1)])))


# ------------------------------------------------------------------ probes and line contraction
def probe_mat(name):
    """symmetric graviton probe sym(nameU, nameW)."""
    return TT.sym(U.V({name + "U": 1}), U.V({name + "W": 1}))


def probe_vec(name):
    return U.V({name: 1})


_SYM = re.compile(r"^(H|T|s)_(.+)$")


def _parse_dot(sym, probes):
    """symbol -> (kind, (a, b)) with a, b basis names; s_x == l.x."""
    m = _SYM.match(sym.name)
    if not m:
        return None
    kind, rest = m.groups()
    if kind == "s":
        a, b = "l", rest
    else:
        a, b = rest.split("_")
    if a in probes or b in probes:
        return kind, a, b
    return None


def contract(expr, pairs):
    """contract probe vectors pairwise with the 4S metric: pairs = [(p, q)] meaning sum_mu p^mu q_mu over a basis.
    Graviton lines use probes (XU, XW) and (YU, YW) with the propagator numerator applied by the caller as pairs.
    Implementation: every monomial is linear in each probe; probe dots are rewritten as index factors and evaluated."""
    probes = {p for pr in pairs for p in pr}
    expr = sp.expand(expr)
    syms = [s for s in expr.free_symbols if (_parse_dot(s, probes) or (None,))[0] == "T"]
    expr = sp.expand(expr.xreplace({s: 0 for s in syms}))   # 4S dots appear as H + T: keep the H copy as the full dot
    lab = {}
    for p, q in pairs:
        lab[p] = lab[q] = L()
    out = 0
    for term in sp.Add.make_args(expr):
        coef = term
        facs = []
        seen = {}
        for s, e in term.as_powers_dict().items():
            pd = _parse_dot(s, probes) if isinstance(s, sp.Symbol) else None
            if pd is None:
                continue
            _, a, b = pd
            for _ in range(int(e)):
                coef = coef / s
                ina, inb = a in probes, b in probes
                if ina and inb:
                    facs.append(("m", ETA(), lab[a], lab[b]))
                elif ina:
                    facs.append(("v", U.V({b: 1}), lab[a]))
                else:
                    facs.append(("v", U.V({a: 1}), lab[b]))
        out += TT.evaluate((coef, facs))
    return sp.expand(out)


def graviton_line(expr, A, B, Z, c=sp.Rational(1, 2)):
    """close a graviton line between probe A (at one vertex) and probe B (other vertex), momentum-independent
    numerator only: [ (eta_ac eta_bd + eta_ad eta_bc)/2 - c eta_ab eta_cd ] / Z ; the i/k^2 is put in by the caller."""
    t1 = contract(expr, [(A + "U", B + "U"), (A + "W", B + "W")])
    t2 = contract(expr, [(A + "U", B + "W"), (A + "W", B + "U")])
    t3 = contract(expr, [(A + "U", A + "W"), (B + "U", B + "W")])
    return sp.expand((sp.Rational(1, 2) * (t1 + t2) - c * t3) / Z)


def ghost_line(expr, A, B, Zg=1):
    return sp.expand(contract(expr, [(A, B)]) / Zg)


# ------------------------------------------------------------------ antifield vertices for s c, s phi (XVI-P P13 s_images)
def V_cs_c_c(zv, c1, q1, c2, q2):
    """i x coefficient of c*(zv) c(c1, q1) c(c2, q2) in int c*_mu (kappa/2) c^r d_r c^mu  (ordering c* c1 c2;
    the c1 <-> c2 exchange carries the Grassmann sign -1)."""
    t = sp.Rational(1, 2) * ((-I) * U.dot(c1, q2) * U.dot(zv, c2) - (-I) * U.dot(c2, q1) * U.dot(zv, c1))
    return sp.expand(I * t)


def V_phis_c_phi(cv, qphi):
    """i x coefficient of phi* c(cv) phi(qphi) in int phi* (kappa/2) c^r d_r phi."""
    return sp.expand(I * sp.Rational(1, 2) * (-I) * U.dot(cv, qphi))


# ------------------------------------------------------------------ fast path: contract two vertices pattern-wise
def _split_probe(expr, probes):
    """expr -> {probe monomial (H-copies only): coefficient}  (T-copies of probe dots dropped, as in contract())."""
    expr = sp.expand(expr)
    drop = [s for s in expr.free_symbols if (_parse_dot(s, probes) or (None,))[0] == "T"]
    expr = sp.expand(expr.xreplace({s: 0 for s in drop}))
    out = {}
    for term in sp.Add.make_args(expr):
        pat, coef = sp.Integer(1), sp.Integer(1)
        for f, e in term.as_powers_dict().items():
            if isinstance(f, sp.Symbol) and _parse_dot(f, probes):
                pat *= f ** e
            else:
                coef *= f ** e
        out[pat] = out.get(pat, 0) + coef
    return out


def contract_vertices(A, B, lines):
    """sum over internal lines between vertex expressions A and B.
    lines: list of ('h', probeA, probeB) graviton (numerator (I - J/2)) or ('c', probeA, probeB, factor) ghost/vector.
    Equivalent to graviton_line/ghost_line applied to sp.expand(A*B), but pattern-wise (memoised)."""
    probes = set()
    for ln in lines:
        if ln[0] == "h":
            probes |= {ln[1] + "U", ln[1] + "W", ln[2] + "U", ln[2] + "W"}
        else:
            probes |= {ln[1], ln[2]}
    SA, SB = _split_probe(A, probes), _split_probe(B, probes)
    cache = {}
    tot = 0
    for pa, ca in SA.items():
        for pb, cb in SB.items():
            key = pa * pb
            if key not in cache:
                v = key
                for ln in lines:
                    v = graviton_line(v, ln[1], ln[2], 1) if ln[0] == "h" else ghost_line(v, ln[1], ln[2], 1 / ln[3])
                cache[key] = v
            if cache[key] != 0:
                tot += cache[key] * ca * cb
    return sp.expand(tot)
