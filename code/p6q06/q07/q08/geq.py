#!/usr/bin/env python3
"""XVI-Q Q08: ghost-equation map on symbol polynomials.  Gamma_{cbar c h} = rho * Gamma_{K c h}[K -> T(cbar)],
T(cbar) = sym(k, cbar) - (1/2)(k.cbar) eta, k = p2 (K / cbar momentum).  Linear in K: X(K) -> bilinear B(u, w)
(rows.polarise), then X(sym(k, gb)) by u = k, w = gb, and X(eta) = sum_mu B(e_mu, e_mu) by index contraction:
H_u_x H_w_y -> H_x_y, T_u_x T_w_y -> T_x_y, mixed -> 0, H_u_w -> 4, T_u_w -> 0 (d -> 4 residue level; hat/tilde
orthogonal).  Validated against the tensor-machinery map of Q05 ghost_eq3 on the loop poles."""
from __future__ import annotations
import sys
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "p6q05"))
import rows as R  # noqa: E402
import uvpole as U  # noqa: E402


def _B(X):
    return sp.expand(R.polarise(sp.expand(X), "aK", "u", "w"))


def _sub_vec(e, name, V):
    rep = {}
    for s in e.free_symbols:
        p = s.name.split("_")
        if len(p) != 3 or name not in p[1:]:
            continue
        o = p[2] if p[1] == name else p[1]
        if o == name:
            continue
        rep[s] = sum(c * (U.H(k, o) if p[0] == "H" else (0 if o.startswith("p") or k.startswith("p") else U.T(k, o)))
                     for k, c in V.items())
    return e.xreplace(rep)


def _trace(B):
    out = 0
    for t in sp.Add.make_args(B):
        c, fu, fw, uw = t, None, None, None
        for f in sp.Mul.make_args(t):
            b, n = f.as_base_exp()
            if not isinstance(b, sp.Symbol):
                continue
            p = b.name.split("_")
            if len(p) != 3:
                continue
            if set(p[1:]) == {"u", "w"}:
                uw = p[0]
                c = c / f
            elif "u" in p[1:]:
                fu = (p[0], p[2] if p[1] == "u" else p[1])
                c = c / f
            elif "w" in p[1:]:
                fw = (p[0], p[2] if p[1] == "w" else p[1])
                c = c / f
        if uw is not None:
            out += c * (4 if uw == "H" else 0)
        elif fu and fw:
            if fu[0] != fw[0]:
                continue
            if fu[0] == "T" and (fu[1].startswith("p") or fw[1].startswith("p")):
                continue
            out += c * (U.H(fu[1], fw[1]) if fu[0] == "H" else U.T(fu[1], fw[1]))
        else:
            raise ValueError(t)
    return sp.expand(out)


def Tmap(X):
    B = _B(X)
    s = _sub_vec(_sub_vec(B, "u", {"p2": sp.Rational(1, 2)}), "w", {"gb2": 1})
    s2 = _sub_vec(_sub_vec(B, "u", {"gb2": sp.Rational(1, 2)}), "w", {"p2": 1})
    return sp.expand(s + s2 - sp.Rational(1, 2) * U.H("gb2", "p2") * _trace(B))


if __name__ == "__main__":
    K = R.strip_t(R.load("G03") + R.load("G07"))
    C = R.strip_t(R.load("G04") + R.load("G09"))
    TK = Tmap(K)
    print("ratio", sp.cancel(C / TK) if TK != 0 else "0", "ratio G04/G03", sp.cancel(R.strip_t(R.load("G04")) / Tmap(R.strip_t(R.load("G03")))))
