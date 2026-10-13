#!/usr/bin/env python3
"""XVI-Q Q07 goal 2: row h.h.c (c(gc1) p1, h(a2) p2, h(a3) p3) of (S, Gamma_1) at kappa^3, UV-pole level.
Channels (linearised antibracket, 'i x coefficient' normalisation as in the Q05/Q06 rows):
  A  Gamma_hhh [h(a1) -> s0 h(c)]                         (G01+G02+G08 + scalar G14+G15, own)
  B  Gamma_hh(h_i, s1 h(h_j, c)),  s1 h = (1/2) L_c h     (own: graviton + ghost + scalar loops)
  C  delta S_EH3/delta h|_{h2 h3} . Gamma_{K c}           (D13 own)
  D  delta S_2/delta h|_{h_i} (EH + gf) . Gamma_{K h_j c} (G03+G07)
  E  s cbar = sigma_F F(h_i), F = D[h]:  F . Gamma_{cbar c h_j}   (G04+G09)  -- channel absent from the Q04 ledger
Tree calibration: (S,S)/2 at kappa^1 on h.h.c = A_t + B_t + sigma_F E_t = 0 fixes sigma_F; C, D enter with the
c.phi.phi K-channel convention (+1).  Polarisation matrices inside loops use t-FREE momenta."""
from __future__ import annotations
import json, sys, time  # noqa: E401
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
Q = HERE.parents[2]
sys.path.insert(0, str(HERE.parents[1] / "p6q05"))
sys.path.insert(0, str(HERE.parent))
import gravity as G  # noqa: E402
import graphs as GR  # noqa: E402
import rows as R  # noqa: E402
import uvpole as U  # noqa: E402
import tensor as TT  # noqa: E402

I = sp.I
p1, p2, l, p3 = GR.p1, GR.p2, GR.l, GR.p3
F1, F2 = U.V({"p1": 1}), U.V({"p2": 1})          # t-free copies (for polarisation-side momenta)
F3 = -(F1 + F2)
LOG = Q / "build" / "logs"
GC = U.vec("gc1")
EXT = {"2": ("a2", p2, F2), "3": ("a3", p3, F3)}


def ext(a):
    return GR.ext(a)


def strip(e):
    return R.strip_t(e)


def subst_vec(e, name, V):
    """linear substitution of the vector polarisation `name` by V = {basis: coef} in H_/T_ symbols."""
    rep = {}
    for s in e.free_symbols:
        parts = s.name.split("_")
        if len(parts) != 3 or parts[0] not in ("H", "T") or name not in parts[1:]:
            continue
        o = parts[2] if parts[1] == name else parts[1]
        if parts[0] == "H":
            rep[s] = sum(c * U.H(k, o) for k, c in V.items())
        else:
            rep[s] = 0 if o in ("p1", "p2") else sum(c * U.T(k, o) for k, c in V.items() if k not in ("p1", "p2"))
    return sp.expand(e.xreplace(rep))


# ------------------------------------------------------------------ tree calibration
def tree():
    Es = TT.sym(p1, GC, -I)
    A = G.V3_grav((Es, p1), (ext("a2"), p2), (ext("a3"), p3))
    B = 0
    E = 0
    for i, j in (("2", "3"), ("3", "2")):
        ai, ki, _ = EXT[i]
        aj, kj, _ = EXT[j]
        M = G.Lc_h(ext(aj), kj, GC, p1)
        M.items = [(sp.Rational(1, 2) * c, u, w) for c, u, w in M.items]
        M.eta = sp.Rational(1, 2) * M.eta
        eh, gf = G.S2_EH_gf(ext(ai), M, ki)
        B += I * (eh + gf)
        E += G.V_cb_c_h(G.Dvec(ext(ai), ki), GC, p1, ext(aj), kj)
    A, B, E = strip(A), strip(B), strip(E)
    sig = sp.Symbol("sigma_F")
    tot = sp.expand(A + B + sig * E)
    gens = sorted([s for s in tot.free_symbols if s != sig], key=str)
    sol = sp.solve(sp.Poly(tot, *gens).coeffs(), sig, dict=True)
    return {"sigma_F": str(sol[0][sig]) if sol else None, "closes": bool(sol),
            "without_E_zero": sp.expand(A + B) == 0, "A_terms": len(sp.Add.make_args(A))}


# ------------------------------------------------------------------ own loops
def gamma_hh(E1, k, E2):
    """own Gamma_hh pole: legs (E1, k) and (E2, -k); graviton bubble (1/2) + ghost loop (-1) + scalar bubble (1/2)."""
    A = G.V3_grav((E1, k), (G.probe_mat("hA1"), -l), (G.probe_mat("hA2"), l - k))
    B = G.V3_grav((E2, -k), (G.probe_mat("hB1"), l), (G.probe_mat("hB2"), k - l))
    n = G.contract_vertices(A, B, [("h", "hA1", "hB1"), ("h", "hA2", "hB2")])
    grav = U.loop_pole(n, [U.V(), -k], 4, sp.Rational(1, 2))
    gA = G.V_cb_c_h(G.probe_vec("gAb"), G.probe_vec("gAc"), l, E1, k)
    gB = G.V_cb_c_h(G.probe_vec("gBb"), G.probe_vec("gBc"), l + k, E2, -k)
    num = G.ghost_line(G.ghost_line(sp.expand(gA * gB), "gAb", "gBc", 1 / GR.GH), "gBb", "gAc", 1 / GR.GH)
    gh = U.loop_pole(-num, [U.V(), k], 4, 1)
    sA = G.V_hphiphi(E1, l, -l - k)
    sB = G.V_hphiphi(E2, l + k, -l)
    sc = U.loop_pole(sp.expand(sA * sB), [U.V(), k], 4, sp.Rational(1, 2))
    return GR.r4(grav, "Q07 hh grav"), GR.r4(gh, "Q07 hh ghost"), GR.r4(sc, "Q07 hh scalar")


def scalar_hhh():
    """own scalar loops for Gamma_hhh: triangle (G14) + V_hh_phi_phi bubbles (G15, 3 channels, symmetry 1/2)."""
    vA = G.V_hphiphi(ext("a1"), l, -l - p1)
    vB = G.V_hphiphi(ext("a2"), l + p1, -l - p1 - p2)
    vC = G.V_hphiphi(ext("a3"), l + p1 + p2, -l)
    tri = U.loop_pole(sp.expand(vA * vB * vC), [U.V(), p1, p1 + p2], 4, 1)
    bub = 0
    legs = {"1": ("a1", p1), "2": ("a2", p2), "3": ("a3", p3)}
    for c in "123":
        a, k = legs[c]
        o = [legs[x] for x in "123" if x != c]
        v1 = G.V_hphiphi(ext(a), l, -l - k)
        v2 = G.V_hhphiphi(ext(o[0][0]), ext(o[1][0]), l + k, -l)
        bub += U.loop_pole(sp.expand(v1 * v2), [U.V(), k], 4, sp.Rational(1, 2))
    return GR.r4(tri, "Q07 G14"), GR.r4(bub, "Q07 G15")


def load_any(g):
    for f in (LOG / f"q06_{g}.json", LOG / f"q05_{g}.json", LOG / f"q07_{g}.json"):
        if f.is_file():
            return sp.sympify(json.loads(f.read_text())["pole_over_U_d4_residue"])
    if g == "G02" and all((LOG / f"q06_G02_{c}.json").is_file() for c in "123"):   # fresh builds: sum of channels
        return sp.expand(sum(load_any(f"G02_{c}") for c in "123"))
    raise FileNotFoundError(g)


def chan_A():
    hhh = sum(load_any(g) for g in ("G01", "G02", "G08", "G14", "G15"))
    return strip(-I * R.polarise(sp.expand(hhh), "a1", "p1", "gc1"))


def chan_B(i):
    j = "3" if i == "2" else "2"
    ai, ki, _ = EXT[i]
    aj, _, Fj = EXT[j]
    M = G.Lc_h(ext(aj), Fj, GC, F1)
    M.items = [(sp.Rational(1, 2) * c, u, w) for c, u, w in M.items]
    M.eta = sp.Rational(1, 2) * M.eta
    g, h, s = gamma_hh(ext(ai), ki, M)
    return {"grav": strip(g), "ghost": strip(h), "scalar": strip(s)}


def chan_C():
    V = G.V3_grav((G.probe_mat("E"), p1), (ext("a2"), p2), (ext("a3"), p3))
    return strip(R.contract_id(R.D13_own("gc1", p1), V, "K", "E", "h"))


def _relabel(e, i):
    """stored X c h poles: c p1, X p2, h(a3) p3.  i = '2': as stored; i = '3': X at p3, h(a2) at p2."""
    return e if i == "2" else R.remap(e, {"p1": R.P1, "p2": R.P3, "a3": "a2"})


def chan_D(i):
    ai, ki, _ = EXT[i]
    Gk = _relabel(load_any("G03") + load_any("G07"), i)
    Gk = R.polarise(sp.expand(Gk), "aK", "KU", "KW")
    eh, gf = G.S2_EH_gf(ext(ai), G.probe_mat("E"), ki)
    return strip(R.contract_id(Gk, I * (eh + gf), "K", "E", "h"))


def chan_E(i):
    ai, _, Fi = EXT[i]
    Gb = _relabel(load_any("G04") + load_any("G09"), i)
    return strip(subst_vec(sp.expand(Gb), "gb2", G.Dvec(ext(ai), Fi)))


TASKS = {"tree": tree, "A": chan_A, "B2": lambda: chan_B("2"), "B3": lambda: chan_B("3"), "C": chan_C,
         "D2": lambda: chan_D("2"), "D3": lambda: chan_D("3"), "E2": lambda: chan_E("2"), "E3": lambda: chan_E("3"),
         "S": scalar_hhh}


if __name__ == "__main__":
    k = sys.argv[1]
    t0 = time.time()
    r = TASKS[k]()
    if k == "S":
        for g, v in zip(("G14", "G15"), r):
            (LOG / f"q07_{g}.json").write_text(json.dumps({"graph": "P6-" + g, "pole_over_U_d4_residue": str(v),
                                                           "n_terms": len(sp.Add.make_args(v))}) + "\n")
        r = {"G14": len(sp.Add.make_args(r[0])), "G15": len(sp.Add.make_args(r[1]))}
    out = {k2: str(v) for k2, v in r.items()} if isinstance(r, dict) else {"value": str(r)}
    (LOG / f"q07_hhc_{k}.json").write_text(json.dumps(out, indent=1) + "\n")
    print(k, {a: (b[:80] if isinstance(b, str) else b) for a, b in out.items()}, round(time.time() - t0, 1))
