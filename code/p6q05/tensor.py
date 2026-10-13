#!/usr/bin/env python3
"""XVI-Q Q05 index-contraction evaluator on top of q04 uvpole vectors.

Factors: ('v', V, i)      vector with index label i
         ('m', Mat, i, j) rank-2 tensor   Mat = eta*eta_4S + sum c u (x) w   (NOT symmetrised; order u -> i, w -> j)
A term is (coef, [factors]); every label appears once (free) or twice (contracted, metric = eta_4S).
Contracted labels resolve into chains (vector ... vector) and cycles (traces); tr eta_4S = 4 (Q4S vertex indices).
Free labels are not supported (vertices are always fully contracted with polarisations / probe vectors)."""
from __future__ import annotations
import sys
from pathlib import Path
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "q04"))
import uvpole as U  # noqa: E402

TR_ETA = 4


class Mat:
    def __init__(self, eta=0, items=()):
        self.eta = eta
        self.items = list(items)

    def T(self):
        return Mat(self.eta, [(c, w, u) for c, u, w in self.items])

    def left(self, v):                    # v . M  (contract v with first index)
        out = U.V()
        if self.eta != 0:
            out = out + self.eta * v
        for c, u, w in self.items:
            out = out + (c * U.dot(v, u)) * w
        return out

    def tr(self):
        return sp.expand(TR_ETA * self.eta + sum(c * U.dot(u, w) for c, u, w in self.items))


def sym(u, w, c=1):
    """c * (u w + w u)/2"""
    h = sp.Rational(1, 2) * c
    return Mat(0, [(h, u, w), (h, w, u)])


def eps_mat(e):
    """q04 Eps (sum c (u w + w u)/2) -> Mat."""
    m = Mat(getattr(e, "eta", 0))
    for c, u, w in e:
        m.items += [(sp.Rational(1, 2) * c, u, w), (sp.Rational(1, 2) * c, w, u)]
    return m


def ETA():
    return Mat(1)


def _chain_value(start_vec, mats):
    v = start_vec
    for M in mats:
        v = M.left(v)
    return v


def evaluate(term):
    coef, facs = term
    # adjacency: label -> list of (factor index, slot)
    occ = {}
    for fi, f in enumerate(facs):
        for s, lab in enumerate(f[2:]):
            occ.setdefault(lab, []).append((fi, s))
    for lab, o in occ.items():
        if len(o) != 2:
            raise ValueError(f"label {lab} appears {len(o)} times")
    used = set()
    total = sp.Integer(1) * coef

    def other_end(fi, s):
        lab = facs[fi][2 + s]
        a, b = occ[lab]
        return b if a == (fi, s) else a

    def oriented(fi, enter_slot):
        """matrix factor fi entered through slot enter_slot -> (Mat oriented, exit slot)."""
        M = facs[fi][1]
        return (M, 1) if enter_slot == 0 else (M.T(), 0)

    # chains from vectors
    for fi, f in enumerate(facs):
        if f[0] != "v" or fi in used:
            continue
        used.add(fi)
        mats = []
        cur = other_end(fi, 0)
        while facs[cur[0]][0] == "m":
            used.add(cur[0])
            M, ex = oriented(cur[0], cur[1])
            mats.append(M)
            cur = other_end(cur[0], ex)
        used.add(cur[0])
        total = total * U.dot(_chain_value(f[1], mats), facs[cur[0]][1])
    # cycles of matrices
    for fi, f in enumerate(facs):
        if fi in used:
            continue
        mats = []
        cur = (fi, 0)
        while True:
            used.add(cur[0])
            M, ex = oriented(cur[0], cur[1])
            mats.append(M)
            nxt = other_end(cur[0], ex)
            if nxt[0] == fi:
                break
            cur = nxt
        total = total * _trace(mats)
    return sp.expand(total)


def _trace(mats):
    if len(mats) == 1:
        return mats[0].tr()
    M0, rest = mats[0], mats[1:]
    out = M0.eta * _trace([Mat(rest[0].eta, rest[0].items)] + rest[1:]) if M0.eta != 0 else 0
    if M0.eta != 0 and len(rest) == 1:
        out = M0.eta * rest[0].tr()
    for c, u, w in M0.items:
        out += c * U.dot(_chain_value(w, rest), u)
    return sp.expand(out)


def ev(*terms):
    return sp.expand(sum(evaluate(t) for t in terms))
