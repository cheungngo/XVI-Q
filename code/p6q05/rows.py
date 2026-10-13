#!/usr/bin/env python3
"""XVI-Q Q05 kappa^3 identity rows at UV-pole level: (S, Gamma_1^pole) restricted to a three-field monomial.
All amplitudes are 'i x coefficient' (vertices) / 'i Gamma' (loop poles, units U) in the same normalisation, so a
row is linear in them.  The relative antibracket sign sigma of each antifield channel is CALIBRATED on the classical
identity (S_0, S_1) = 0 (tree level, kappa^1) -- the O.2-verified master equation -- and then used unchanged at one loop.
Row c.phi.phi (c p1, phi p2, phi p3):
   (S0, G3): Gamma_h phi phi [h -> s0 h(c)]  +  sigma * sum_phi p_phi^2 Gamma_{phi* c phi}
   (S1, G2): delta S_1/delta h |phi phi  contracted with Gamma_{K c} (D13 own pole)      (Gamma_phi phi = 0, Q02)
Pole-level requirement: every evanescent structure dies under the s0 polarisation (T(p1, c) = 0), so the row must
vanish EXACTLY (no evanescent counterterm can absorb a remainder)."""
from __future__ import annotations
import json, sys  # noqa: E401
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
Q = HERE.parents[1]
sys.path.insert(0, str(HERE))
import gravity as G  # noqa: E402
import uvpole as U  # noqa: E402
import tensor as TT  # noqa: E402
import graphs as GR  # noqa: E402
import dred  # noqa: E402

I = sp.I
p1, p2, l = GR.p1, GR.p2, GR.l
p3 = GR.p3
LOG = Q / "build" / "logs"


def load(g):
    return sp.sympify(json.loads((LOG / f"q05_{g}.json").read_text())["pole_over_U_d4_residue"])


def polarise(Qa, a, u, w):
    """quadratic form Q(a) in symbols H_a_*, T_a_* -> bilinear B(u, w) (= value at E = sym(u, w))."""
    rep = {}
    for s in Qa.free_symbols:
        n = s.name
        if not (n.startswith("H_") or n.startswith("T_")):
            continue
        kind, x, y = n.split("_")
        if a not in (x, y):
            continue
        F = U.H if kind == "H" else U.T
        if x == a and y == a:
            rep[s] = sum(F(v1, v2) for v1 in ("tu", "tw") for v2 in ("tu", "tw"))
        else:
            o = y if x == a else x
            rep[s] = F("tu", o) + F("tw", o)
    e = sp.expand(Qa.xreplace(rep))
    tu = [s for s in e.free_symbols if "tu" in s.name.split("_")]
    tw = [s for s in e.free_symbols if "tw" in s.name.split("_")]
    mixed = 0
    for term in sp.Add.make_args(e):
        nu = sum(sp.degree(term, s) for s in tu)
        nw = sum(sp.degree(term, s) for s in tw)
        if nu == 1 and nw == 1:
            mixed += term
    B = sp.expand(mixed / 2)
    # substitute u, w
    rep2 = {}
    for s in B.free_symbols:
        kind, x, y = (s.name.split("_") + [None])[:3] if s.name.count("_") == 2 else (None, None, None)
        if kind is None or ("tu" not in (x, y) and "tw" not in (x, y)):
            continue
        sub = {"tu": u, "tw": w}
        X = sub.get(x, x)
        Y = sub.get(y, y)
        rep2[s] = _dot_sym(kind, X, Y)
    return sp.expand(B.xreplace(rep2))


def _dot_sym(kind, X, Y):
    mom = lambda z: z.startswith("p")  # noqa: E731
    if kind == "T" and (mom(X) or mom(Y)):
        return 0
    if X == Y and mom(X) and X == "p1":
        return U.H("p1", "p1")
    return (U.H if kind == "H" else U.T)(X, Y)


def kc_tree(KE, gc, q):
    """tree K.c amplitude: i (1/2) K:(d c + d c),  d -> -i q."""
    M = TT.sym(q, gc, -I)
    i1, i2 = G.L(), G.L()
    return sp.expand(I * TT.ev((1, [("m", KE, i1, i2), ("m", M, i2, i1)])))


def k_channel(Kamp, phis):
    """delta S1/delta h|phi phi (vertex V_h_phi_phi) contracted with (K-amplitude / i) through the identity pairing."""
    V = G.V_hphiphi(G.probe_mat("E"), *phis)
    e = sp.expand(Kamp * V / I)
    t1 = G.contract(e, [("KU", "EU"), ("KW", "EW")])
    t2 = G.contract(e, [("KU", "EW"), ("KW", "EU")])
    return sp.expand((t1 + t2) / 2)


def phistar_tree(gc):
    return sp.expand(U.dot(p3, p3) * G.V_phis_c_phi(gc, p2) + U.dot(p2, p2) * G.V_phis_c_phi(gc, p3))


def calibrate_cphiphi():
    gc = U.vec("gc1")
    KE = G.probe_mat("K")
    a = k_channel(kc_tree(KE, gc, p1), (p2, p3))
    b = phistar_tree(gc)
    r = sp.cancel(-a / b)
    return ([r] if not r.free_symbols else []), sp.expand(a), sp.expand(b)


def D13_own_K():
    """own D13 K.c pole as a function of a probe K (KU, KW), c pol gc1 momentum p1 (K carries -p1)."""
    gc = U.vec("gc1")
    A = G.V_cb_c_h(G.probe_vec("gA"), gc, p1, G.probe_mat("hA"), l)
    B = G.V_K_c_h(G.probe_mat("K"), G.probe_vec("gB"), p1 + l, G.probe_mat("hB"), -l)
    n = G.contract_vertices(A, B, [("h", "hA", "hB"), ("c", "gA", "gB", GR.GH)])
    return sp.expand(dred.substitute_d(U.loop_pole(n, [U.V(), p1], 3, 1), 4, "Q05 row: D13 own pole"))


def row_cphiphi(sig):
    hpp = sum(load(g) for g in ("G05", "G06", "G12", "G13"))
    t_h = -I * polarise(hpp, "a1", "p1", "gc1")                     # E = s0 h = sym(p1, gc, -i)
    g11a = load("G11")                                               # phi p2, phi* p3
    g11b = GR.G11_swapped()                                          # phi p3, phi* p2
    t_s = sig * (U.dot(p3, p3) * g11a + U.dot(p2, p2) * g11b)
    t_k = k_channel(D13_own_K(), (p2, p3))
    t = sp.Symbol("t")
    t_h, t_s, t_k = (sp.expand(x.subs(t, 1)) for x in (t_h, t_s, t_k))   # t is a grading parameter only
    R = sp.expand(t_h + t_s + t_k)
    return {"S0_hphiphi": str(t_h), "S0_phistar": str(sp.expand(t_s)), "S1_Kc": str(t_k), "residual": str(R),
            "zero": R == 0, "without_S1_Kc_zero": sp.expand(t_h + t_s) == 0}


if __name__ == "__main__":
    sol, a, b = calibrate_cphiphi()
    print("calibration sigma:", sol, "| tree K-channel:", a, "| tree phi*:", b)
    if sol:
        r = row_cphiphi(sol[0])
        print({k: (v if k in ("zero", "without_S1_Kc_zero") else v[:200]) for k, v in r.items()})


# ======================================================================= ghost rows (two identical Grassmann c's)
def _V(name):
    return U.V({name: 1})


def remap(expr, mp):
    """momentum relabelling of a pole expression: mp maps 'p1','p2' to U.V vectors (no t); p3 implicit;
    polarisation renames via mp entries for non-momentum names (str -> str)."""
    rep = {}
    for s in expr.free_symbols:
        parts = s.name.split("_")
        if len(parts) != 3 or parts[0] not in ("H", "T"):
            continue
        kind, x, y = parts
        def vec(z):  # noqa: E306
            m = mp.get(z, z)
            return m if isinstance(m, U.V) else _V(m)
        X, Y = vec(x), vec(y)
        if kind == "T":
            mom = any(k.startswith("p") for k in list(X) + list(Y))
            rep[s] = 0 if mom else U.T(next(iter(X)), next(iter(Y)))
        else:
            rep[s] = sp.expand(sum(cx * cy * U.H(*sorted((kx, ky))) if True else 0
                                   for kx, cx in X.items() for ky, cy in Y.items()))
    return sp.expand(expr.xreplace(rep))


P1, P2 = _V("p1"), _V("p2")
P3 = -(P1 + P2)


def strip_t(e):
    return sp.expand(sp.sympify(e).subs(sp.Symbol("t"), 1))


def contract_id(A, B, pa, pb, kind):
    """identity pairing of probe pa (in A) with probe pb (in B): kind 'h' (sym pair) or 'c' (vector)."""
    e = sp.expand(A * B / I)
    if kind == "c":
        return G.contract(e, [(pa, pb)])
    t1 = G.contract(e, [(pa + "U", pb + "U"), (pa + "W", pb + "W")])
    t2 = G.contract(e, [(pa + "U", pb + "W"), (pa + "W", pb + "U")])
    return sp.expand((t1 + t2) / 2)


def kin_cbc(gb, X, q):
    """tree cbar(gb) c(X, q) amplitude: i * coefficient of i cbar.Box c (Box -> -q^2)."""
    return sp.expand(I * I * (-U.dot(q, q)) * U.dot(gb, X))


def kin_Kc(KE, X, q):
    return kc_tree(KE, X, q)


def D13_own(gcname, q):
    gc = U.vec(gcname)
    A = G.V_cb_c_h(G.probe_vec("gA"), gc, q, G.probe_mat("hA"), l)
    B = G.V_K_c_h(G.probe_mat("K"), G.probe_vec("gB"), q + l, G.probe_mat("hB"), -l)
    n = G.contract_vertices(A, B, [("h", "hA", "hB"), ("c", "gA", "gB", GR.GH)])
    return strip_t(dred.substitute_d(U.loop_pole(n, [U.V(), q], 3, 1), 4, "Q05 row: D13 own pole"))


def D06_own(gbname, X, q):
    """own D06 cbar(gb) c(probe X, momentum q) pole."""
    A = G.V_cb_c_h(G.probe_vec("gA"), G.probe_vec(X), q, G.probe_mat("hA"), l)
    B = G.V_cb_c_h(U.vec(gbname), G.probe_vec("gB"), q + l, G.probe_mat("hB"), -l)
    n = G.contract_vertices(A, B, [("h", "hA", "hB"), ("c", "gA", "gB", GR.GH)])
    return strip_t(dred.substitute_d(U.loop_pole(n, [U.V(), q], 4, 1), 4, "Q05 row: D06 own pole"))


def ghost_row(target, loop=True, sig_c=None, drop=()):
    """target 'ccb' (c p1, c p2, cbar p3) or 'ccK' (c p1, c p2, K p3).  Channels:
       h:  Gamma_{X c h}[h -> s0 h(c_j)]                (tree: S1 vertex; loop: G04+G09 / G03+G07)
       K:  delta S1/delta h |_{X c} . Gamma_{K c}       (loop only; D13 own)
       c:  delta Gamma/delta c . delta Gamma/delta c*   (tree: kinetic x V_cs; loop: kinetic x G10 + Gamma_{X c} x V_cs)
    antisymmetrised in (c1, p1) <-> (c2, p2) with the Grassmann sign -1."""
    t = sp.Symbol("t")
    gc1, gc2 = U.vec("gc1"), U.vec("gc2")
    if target == "ccb":
        Xside = lambda cv, q, E, k: G.V_cb_c_h(U.vec("gb3"), cv, q, E, k)      # noqa: E731
        kin = lambda Xp, q: kin_cbc(U.vec("gb3"), Xp, q)                       # noqa: E731
        loopX = ("G04", "G09")      # computed with c p1, cbar p2, h p3
    else:
        Xside = lambda cv, q, E, k: G.V_K_c_h(GR.ext("aK"), cv, q, E, k)        # noqa: E731
        kin = lambda Xp, q: kin_Kc(GR.ext("aK"), Xp, q)                         # noqa: E731
        loopX = ("G03", "G07")
    p1v, p2v = U.vec("p1"), U.vec("p2")
    ch = {"h": 0, "K": 0, "c": 0}
    for (ci, pi, cj, pj, sgn) in (("gc1", p1v, "gc2", p2v, 1), ("gc2", p2v, "gc1", p1v, -1)):
        # h channel: X(c_i, p_i) with h -> s0 h(c_j, p_j)
        Es = TT.sym(pj, U.vec(cj), -I)
        if not loop:
            ch["h"] += sgn * Xside(U.vec(ci), pi, Es, pj)
        else:
            # loop Gamma_{X c h} stored with c(gc1,p1), X at p2, h(a3) at p3: relabel (p1,p2,p3) -> (p_i, p3, p_j)
            Gx = sum(load(g) for g in loopX)
            Pi, Pj = (P1, P2) if ci == "gc1" else (P2, P1)
            Gx = remap(Gx, {"p1": Pi, "p2": P3, "gc1": ci, "gb2": "gb3"})
            ch["h"] += sgn * (-I) * polarise(Gx, "a3", "p" + ("1" if pj is p1v else "2"), cj)
            # K channel: delta S1/delta h |_{X c_j} (h momentum p_i) . Gamma_{K c}(c_i, p_i)
            Kamp = D13_own(ci, pi)
            ch["K"] += sgn * contract_id(Kamp, Xside(U.vec(cj), pj, G.probe_mat("E"), pi), "K", "E", "h")
    # c channel (already antisymmetric through V_cs_c_c)
    Vcs = G.V_cs_c_c(G.probe_vec("Z"), gc1, p1v, gc2, p2v)
    q3 = -(p1v + p2v)   # momentum of the c leg contracted into X: c(X, q) with q = -p3 = p1 + p2
    ch["c"] += contract_id(kin(G.probe_vec("X"), q3), Vcs, "X", "Z", "c")
    if loop:
        ch["c"] = 0                                                                   # the tree piece is not a loop term
        g10 = load("G10")                                                            # c* pol z at p3
        kz = kin(G.probe_vec("X"), q3)
        # Gamma_{c* c c}(z) . kinetic(X): rename the c* polarisation z -> probe ZZ and pair it with X
        g10p = g10.xreplace({x: sp.Symbol("_".join("ZZ" if q == "z" else q for q in x.name.split("_")))
                             for x in g10.free_symbols if "z" in x.name.split("_")})
        ch["c1"] = sp.expand(G.contract(sp.expand(kz * g10p / I), [("X", "ZZ")]))
        if target == "ccb":
            gx = D06_own("gb3", "X", q3)
        else:
            gx = strip_t(D13_K_ext("X", q3))
        ch["c2"] = contract_id(gx, Vcs, "X", "Z", "c")
    ch = {k: strip_t(v) for k, v in ch.items()}
    if sig_c is None:
        return ch
    R = sp.expand(sum(v * (sig_c if k in ("c", "c1", "c2") else 1) for k, v in ch.items() if k not in drop))
    return ch, R


def D13_K_ext(Xc, q):
    """own D13 with the EXTERNAL K = ext('aK') and the c leg a probe Xc (momentum q)."""
    A = G.V_cb_c_h(G.probe_vec("gA"), G.probe_vec(Xc), q, G.probe_mat("hA"), l)
    B = G.V_K_c_h(GR.ext("aK"), G.probe_vec("gB"), q + l, G.probe_mat("hB"), -l)
    n = G.contract_vertices(A, B, [("h", "hA", "hB"), ("c", "gA", "gB", GR.GH)])
    return dred.substitute_d(U.loop_pole(n, [U.V(), q], 3, 1), 4, "Q05 row: D13 own pole (ext K)")
