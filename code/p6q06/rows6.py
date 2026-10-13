#!/usr/bin/env python3
"""XVI-Q Q06 ghost rows c.c.cbar / c.c.K at UV-pole level -- diagnosis and fix of the Q05 OPEN residual.
Channels per orientation o in {1: (c_i, c_j) = (c1, c2), 2: (c2, c1)} as in p6q05/rows.ghost_row:
  h_o : Gamma_{X c_i h}[h -> s0 h(c_j)]         (G03+G07 / G04+G09)
  K_o : delta S1/delta h|_{X c_j} . Gamma_{K c_i} (D13 own)
  c1  : Gamma_{c* c c} (G10) x kinetic X.c ;  c2 : Gamma_{X c} (D13 / D06 own) x V_{c* c c}
FIX (Q06): in the K channel the external c_i (inside Gamma_{K c_i}) stands to the RIGHT of c_j (inside the S1 vertex),
the opposite order to the h channel (c_i in Gamma, c_j in s0 h); bringing the monomial to the row ordering c1 c2 costs
one extra Grassmann transposition: K_o enters with sign -sgn(o) where Q05 used +sgn(o).  The tree row has no K channel,
so the tree calibration (sigma_c) is unchanged.
Result: with the fix the 4D (non-evanescent) part of both rows vanishes; the channel coefficients are the UNIQUE
solution of an unconstrained 6-coefficient fit (each channel and orientation free), so every channel is needed.
What remains is a purely evanescent remainder (every term carries one T(pol, pol)), recorded explicitly."""
from __future__ import annotations
import hashlib, json, sys  # noqa: E401
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
Q = HERE.parents[1]
sys.path.insert(0, str(HERE.parent / "p6q05"))
import rows as R  # noqa: E402
import gravity as G  # noqa: E402
import uvpole as U  # noqa: E402
import graphs as GR  # noqa: E402

I = sp.I
SIGMA_C = {"ccb": -1, "ccK": 1}          # tree calibration (Q05, unchanged)


def parts(target):
    gc1, gc2 = U.vec("gc1"), U.vec("gc2")
    p1v, p2v = U.vec("p1"), U.vec("p2")
    if target == "ccb":
        def Xside(cv, q, E, k):
            return G.V_cb_c_h(U.vec("gb3"), cv, q, E, k)

        def kin(Xp, q):
            return R.kin_cbc(U.vec("gb3"), Xp, q)
        loopX = ("G04", "G09")
    else:
        def Xside(cv, q, E, k):
            return G.V_K_c_h(GR.ext("aK"), cv, q, E, k)

        def kin(Xp, q):
            return R.kin_Kc(GR.ext("aK"), Xp, q)
        loopX = ("G03", "G07")
    out = {}
    for (ci, pi, cj, pj, o) in (("gc1", p1v, "gc2", p2v, 1), ("gc2", p2v, "gc1", p1v, 2)):
        Gx = sum(R.load(g) for g in loopX)
        Pi = R.P1 if ci == "gc1" else R.P2
        Gx = R.remap(Gx, {"p1": Pi, "p2": R.P3, "gc1": ci, "gb2": "gb3"})
        out[f"h{o}"] = (-I) * R.polarise(Gx, "a3", "p1" if o == 2 else "p2", cj)
        out[f"K{o}"] = R.contract_id(R.D13_own(ci, pi), Xside(U.vec(cj), pj, G.probe_mat("E"), pi), "K", "E", "h")
    Vcs = G.V_cs_c_c(G.probe_vec("Z"), gc1, p1v, gc2, p2v)
    q3 = -(p1v + p2v)
    g10 = R.load("G10")
    kz = kin(G.probe_vec("X"), q3)
    g10p = g10.xreplace({x: sp.Symbol("_".join("ZZ" if q == "z" else q for q in x.name.split("_")))
                         for x in g10.free_symbols if "z" in x.name.split("_")})
    out["c1"] = sp.expand(G.contract(sp.expand(kz * g10p / I), [("X", "ZZ")]))
    gx = R.D06_own("gb3", "X", q3) if target == "ccb" else R.strip_t(R.D13_K_ext("X", q3))
    out["c2"] = R.contract_id(gx, Vcs, "X", "Z", "c")
    return {k: R.strip_t(v) for k, v in out.items()}


def coeffs(target, q05=False):
    s = SIGMA_C[target]
    kq = 1 if q05 else -1
    return {"h1": 1, "h2": -1, "K1": kq, "K2": -kq, "c1": s, "c2": s}


def four_d(e):
    return sp.expand(e.xreplace({x: 0 for x in e.free_symbols if x.name.startswith("T_")}))


def fit(P):
    cs = sp.symbols("x0:6")
    ks = sorted(P)
    tot = sp.expand(sum(c * four_d(P[k]) for c, k in zip(cs, ks)))
    gens = [x for x in tot.free_symbols if x not in cs]
    sol = sp.linsolve(sp.Poly(tot, *gens).coeffs(), cs)
    return ks, [str(v) for v in list(sol)[0]] if sol else None


def row(target):
    P = parts(target)
    C, Cq = coeffs(target), coeffs(target, q05=True)
    Rf = sp.expand(sum(C[k] * P[k] for k in P))
    Rq = sp.expand(sum(Cq[k] * P[k] for k in P))
    muts = {}
    for k in P:
        muts[f"drop_{k}"] = four_d(sp.expand(Rf - C[k] * P[k])) != 0
    muts["flip_sigma_c"] = four_d(sp.expand(Rf - 2 * C["c1"] * (P["c1"] + P["c2"]))) != 0
    muts["Q05_K_sign"] = four_d(Rq) != 0
    ev = Rf
    ks, sol = fit(P)
    return {"row": target, "channel_coefficients": {k: str(v) for k, v in C.items()},
            "channel_terms": {k: len(sp.Add.make_args(v)) for k, v in P.items()},
            "unconstrained_fit_4D": {"channels": ks, "solution_up_to_scale": sol},
            "residual_4D": str(four_d(Rf)), "residual_4D_zero": four_d(Rf) == 0,
            "evanescent_remainder": str(ev), "evanescent_remainder_n_terms": len(sp.Add.make_args(ev)),
            "evanescent_remainder_sha256": hashlib.sha256(str(ev).encode()).hexdigest(),
            "evanescent_remainder_every_term_has_T": all(any(x.name.startswith("T_") for x in t.free_symbols)
                                                          for t in sp.Add.make_args(ev)),
            "Q05_assembly_residual_4D_n_terms": len(sp.Add.make_args(four_d(Rq))),
            "mutations_detected": muts}


if __name__ == "__main__":
    out = {t: row(t) for t in ("ccK", "ccb")}
    (Q / "build" / "logs" / "q06_ghost_rows.json").write_text(json.dumps(out, indent=1) + "\n")
    for t, r in out.items():
        print(t, r["residual_4D_zero"], r["evanescent_remainder_n_terms"], r["unconstrained_fit_4D"], r["mutations_detected"])
