"""XVI-Q Q01 O.2(b): explicit O(kappa^n) antibracket expansion, af-0 gravitational sector, n = 0, 1, 2, on an
INDEPENDENT jet engine (this file; exact sparse polynomial ring over QQ), from the ell_2..ell_4 exports.

Derivation (HT antibracket, XVI-O convention sheet): S = S_EH + S_af + ..., S_af = int h*^{mn} (s h)_{mn} + ...,
S_EH = sum_n kappa^{n-2} int ell_n,  s = s0 + kappa s1 (+ kappa^2 s2, s2 = 0 DERIVED: P13_brst_derivation).
(S,S)|_{Phi*=0, grav} = 2 int sum_{m<=n} (delta S_EH / delta h_{mn}) (s h)_{mn}.  The kappa^n coefficient
(n = 0,1,2) is 2 int W_n with W_n = sum_{m<=n'} [E_{mn'}(ell_{n+2}) (s0 h)_{mn'} + E_{mn'}(ell_{n+1}) (s1 h)_{mn'}];
for n = 2 this IS the af-0 grav part of 2(S0,S2) + (S1,S1) written in components.  W_n is linear in the ghost, so
int W_n = 0 for all compactly supported c  <=>  E_{c^mu}(W_n) = 0 for mu = 0..3 (exact, no exactness theorem needed:
int c^mu N_mu = 0 for all c <=> N_mu = 0).  Conventions: s0 h_mn = a(eta_nn d_m c^n + eta_mm d_n c^m),
s1 h_mn = a(c^r d_r h_mn + h_rn d_m c^r + h_mr d_n c^r), a = 1/2 (XVI-O sheet: xi = a c, g = eta + kappa h)."""
from __future__ import annotations
import itertools
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.rings import ring

ETA = [1, -1, -1, -1]
A = QQ(1, 2)
IDX = [(a, b) for a in range(4) for b in range(a, 4)]
HJ, CJ = 3, 2                                     # jet orders needed (see docstring of main)


def multis(k):
    return [tuple(sorted(t)) for t in itertools.combinations_with_replacement(range(4), k)]


ALL_H = [(ab, al) for ab in IDX for k in range(HJ + 1) for al in multis(k)]
ALL_C = [(mu, al) for mu in range(4) for k in range(CJ + 1) for al in multis(k)]
NAMES = [f"h{a}{b}_" + "".join(map(str, al)) for (a, b), al in ALL_H] + [f"c{m}_" + "".join(map(str, al)) for m, al in ALL_C]
R, *GENS = ring(",".join(NAMES), QQ)
G = dict(zip(NAMES, GENS))
POS = {g: i for i, g in enumerate(GENS)}


def h(a, b, al=()):
    a, b = sorted((a, b))
    return G[f"h{a}{b}_" + "".join(map(str, tuple(sorted(al))))]


def c(m, al=()):
    return G[f"c{m}_" + "".join(map(str, tuple(sorted(al))))]


SHIFT = {}
for i, nm in enumerate(NAMES):
    base, al = nm.split("_")
    al = tuple(int(x) for x in al)
    for r in range(4):
        nn = base + "_" + "".join(map(str, tuple(sorted(al + (r,)))))
        SHIFT[(i, r)] = G.get(nn)


def D(p, r):
    """total derivative d/dx^r on the jet ring (exact)."""
    out = R.zero
    used = set()
    for mon in p.monoms():
        used.update(i for i, e in enumerate(mon) if e)
    for i in used:
        t = SHIFT[(i, r)]
        if t is None:
            raise RuntimeError(f"jet order exceeded at {NAMES[i]}")
        out += p.diff(GENS[i]) * t
    return out


def euler(p, var):
    """E_var(p) = sum_al (-D)^al dp/d var_al  for var = ('h',a,b) or ('c',m)."""
    out = R.zero
    kind = var[0]
    for k in range(HJ + 1):
        for al in multis(k):
            g = h(var[1], var[2], al) if kind == "h" else (c(var[1], al) if k <= CJ else None)
            if g is None:
                continue
            q = p.diff(g)
            if not q:
                continue
            for r in al:
                q = -D(q, r)
            out += q
    return out


def load_export(path):
    e = sp.sympify(Path(path).read_text())
    out = R.zero
    for term in sp.Add.make_args(sp.expand(e)):
        coeff, rest = term.as_coeff_Mul()
        m = R(QQ(Fraction(str(coeff))))
        for f in sp.Mul.make_args(rest):
            base, ex = f.as_base_exp()
            n = base.name
            if n.startswith("dH"):
                g = h(int(n[2]), int(n[3]), (int(n.split("_")[1]),))
            else:
                g = h(int(n[1]), int(n[2]))
            m *= g ** int(ex)
        out += m
    return out


def s0(m, n):
    return A * (ETA[n] * c(n, (m,)) + ETA[m] * c(m, (n,)))


def s1(m, n):
    return A * (sum(c(r) * h(m, n, (r,)) for r in range(4)) + sum(h(r, n) * c(r, (m,)) for r in range(4))
                + sum(h(m, r) * c(r, (n,)) for r in range(4)))


def main(exports_dir, mutate=None):
    t0 = time.time()
    ell = {k: load_export(Path(exports_dir) / f"ell{k}_generating_polynomial.srepr") for k in (2, 3, 4)}
    if mutate == "scale_ell4_term":
        t = ell[4].terms()[0]
        ell[4] += R({t[0]: t[1]})
    if mutate == "drop_s1":
        global s1
        s1 = lambda m, n: R.zero  # noqa: E731
    E = {k: {ab: euler(ell[k], ("h",) + ab) for ab in IDX} for k in (2, 3, 4)}
    out = {}
    for n in (0, 1, 2):
        W = sum((E[n + 2][ab] * s0(*ab) for ab in IDX), R.zero)
        if n >= 1:
            W += sum((E[n + 1][ab] * s1(*ab) for ab in IDX), R.zero)
        N = {mu: euler(W, ("c", mu)) for mu in range(4)}
        out[f"kappa{n}"] = {"W_terms": len(W.terms()), "E_c_W_terms": {mu: len(N[mu].terms()) for mu in range(4)},
                            "antibracket_af0_vanishes": all(not N[mu] for mu in range(4))}
    out["seconds"] = round(time.time() - t0, 1)
    return out


if __name__ == "__main__":
    res = {"route": "independent jet engine (XVI-Q); W_n = E(ell_{n+2}) s0 h + E(ell_{n+1}) s1 h; E_c(W_n) = 0",
           "unmutated": main(sys.argv[1])}
    for m in ("scale_ell4_term", "drop_s1"):
        r = main(sys.argv[1], mutate=m)
        res[f"mutation_{m}"] = {"kappa2_vanishes": r["kappa2"]["antibracket_af0_vanishes"],
                                "detected": not r["kappa2"]["antibracket_af0_vanishes"]}
    print(json.dumps(res, indent=1))
