#!/usr/bin/env python3
"""XVI-Q Q06 fast path for the large hhh graphs (G01 V3^3 triangle, G02 V3-V4 bubble): the SAME algorithm as
p6q05/gravity.contract_vertices + q04/uvpole.loop_pole, but the heavy polynomial products are done in a sparse
polynomial ring over QQ(i) (sympy PolyRing) instead of Expr expansion, and every partial product is truncated in the
Taylor parameter t at the pole order omega (exact: t-degrees only add).  d never enters the ring: it appears only in
uvpole.integrate_pole, which is applied to the final t^omega coefficient as an Expr (guarded, as before).
Equivalence with the Expr path is checked on G05/G06/G12 in q06_certs (bit-identical poles)."""
from __future__ import annotations
import itertools, sys  # noqa: E401
from pathlib import Path
import sympy as sp
from sympy.polys.rings import ring
from sympy.polys.domains import QQ_I

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "p6q05"))
import gravity as G  # noqa: E402
import uvpole as U  # noqa: E402

t, Lsym = U.t, U.L
Linv = sp.Symbol("Linv")


class Fast:
    def __init__(self, names):
        names = sorted(set(names))
        gens = [U.H(x, y) for x, y in itertools.combinations_with_replacement(names, 2)]
        gens += [U.T(x, y) for x, y in itertools.combinations_with_replacement(names, 2)]
        gens += [U.S(x) for x in names] + [t, Lsym, Linv]
        self.gens = list(dict.fromkeys(gens))
        self.R, *_ = ring(self.gens, QQ_I)
        self.idx = {g: i for i, g in enumerate(self.gens)}
        self.it = self.idx[t]

    def to(self, e):
        e = sp.expand(sp.sympify(e).xreplace({Lsym ** -k: Linv ** k for k in range(1, 12)}))
        return self.R.from_expr(e) if e != 0 else self.R.zero

    def back(self, p):
        e = p.as_expr()
        return sp.expand(e.xreplace({Linv: 1 / Lsym}))

    def trunc(self, p, w):
        return self.R({m: c for m, c in p.terms() if m[self.it] <= w})

    def split(self, p, probes):
        """{pattern exponent tuple over probe gens: coefficient ring element}; T-copies of probe dots dropped."""
        pg = [i for i, g in enumerate(self.gens) if g.name.startswith(("H_", "s_")) and
              any(q in probes for q in g.name.split("_")[1:])]
        tg = {i for i, g in enumerate(self.gens) if g.name.startswith("T_") and
              any(q in probes for q in g.name.split("_")[1:])}
        out = {}
        for m, c in p.terms():
            if any(m[i] for i in tg):
                continue
            key = tuple((i, m[i]) for i in pg if m[i])
            mm = list(m)
            for i, _ in key:
                mm[i] = 0
            out.setdefault(key, []).append((tuple(mm), c))
        return {k: self.R(dict(v)) for k, v in out.items()}

    def pattern_expr(self, key):
        e = sp.Integer(1)
        for i, k in key:
            e *= self.gens[i] ** k
        return e

    def contract(self, A, B, lines, w):
        probes = set()
        for ln in lines:
            probes |= ({ln[1] + "U", ln[1] + "W", ln[2] + "U", ln[2] + "W"} if ln[0] == "h" else {ln[1], ln[2]})
        SA, SB = self.split(A, probes), self.split(B, probes)
        cache = {}
        tot = self.R.zero
        for ka, ca in SA.items():
            for kb, cb in SB.items():
                key = (ka, kb)
                if key not in cache:
                    v = self.pattern_expr(ka) * self.pattern_expr(kb)
                    for ln in lines:
                        v = G.graviton_line(v, ln[1], ln[2], 1) if ln[0] == "h" else G.ghost_line(v, ln[1], ln[2], 1 / ln[3])
                    cache[key] = self.to(v)
                if cache[key]:
                    tot += self.trunc(self.trunc(cache[key] * ca, w) * cb, w)
        return tot

    def loop_pole(self, numer, dens, omega, prefactor):
        den = self.to(U.taylor_denominators(dens, omega))
        full = self.trunc(numer, omega) * den
        top = self.R({m: c for m, c in full.terms() if m[self.it] == omega})
        e = self.back(top).subs(t, 1) if top else sp.Integer(0)
        return sp.expand(prefactor * sp.I ** len(dens) * U.integrate_pole(e))
