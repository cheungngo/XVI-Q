#!/usr/bin/env python3
"""XVI-Q Q05: validate own gravity/ghost rules against XVI-P pinned two-point poles.
  D03 (ghost loop, hh)  D06 (cbar c: ghost-graviton bubble)  D13 (K'.c)  D01 (graviton loop, V3 x V3, symmetry 1/2)
usage: validate.py NAME  -> build/logs/q05_val_NAME.json  (ratio mine / XVI-P per structure; constant ratio required)"""
from __future__ import annotations
import json, sys, time  # noqa: E401
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
Q = HERE.parents[1]
sys.path.insert(0, str(HERE))
import gravity as G  # noqa: E402
import uvpole as U  # noqa: E402
import tensor as TT  # noqa: E402
import dred  # noqa: E402

PV = Q / "provenance" / "XVI-P" / "certificates"
GHOST_NUM = 1 / (-sp.I)          # <c cbar> = -1/k^2 = (i/k^2) * i  (loop_pole supplies i/k^2)


def res4(e, why):
    return sp.expand(dred.substitute_d(e, 4, why))


def hh_reference(result):
    """XVI-P hh structures contracted with a1 a1 a2 a2, orbit-weighted (the mapping validated on D11 in Q04)."""
    H, T = U.H, U.T
    idx = {"m": 1, "n": 1, "r": 2, "s": 2}
    tot = 0
    for k, v in result.items():
        parts = k.split(" ")
        e = int(parts[-1][-1]) if parts[-1].startswith("(p.p)") else 0
        parts = [f for f in parts if not f.startswith("(p.p)")]
        val = H("p1", "p1") ** e
        orbit = set()
        for s1 in (0, 1):
            for s2 in (0, 1):
                mp = {"m": "nm"[s1], "n": "mn"[s1], "r": "sr"[s2], "s": "rs"[s2]}
                ren = []
                for f in parts:
                    nm, ix = f.split("_")
                    nix = "".join(mp[c] for c in ix)
                    ren.append(nm + "_" + ("".join(sorted(nix)) if nm in ("g4", "gt") else nix))
                orbit.add(tuple(sorted(ren)))
        for f in parts:
            nm, ix = f.split("_")
            a, b = (f"a{idx[ix[0]]}", f"a{idx[ix[1]]}") if len(ix) == 2 else (f"a{idx[ix]}", None)
            val *= {"g4": lambda: H(a, b) + T(a, b), "gt": lambda: T(a, b), "p": lambda: H(a, "p1")}[nm]()
        tot += len(orbit) * sp.sympify(v["UV_pole"]) * val
    return sp.expand(tot)


def ratio(mine, ref):
    r = sp.cancel(mine / ref) if ref != 0 else None
    return {"ratio": str(r), "constant": r is not None and not r.free_symbols}


def ext_mat(a):
    return TT.sym(U.vec(a), U.vec(a))


def D03():
    p, l = U.vec("p1"), U.vec("l")
    E1, E2 = ext_mat("a1"), ext_mat("a2")
    A = G.V_cb_c_h(G.probe_vec("gAb"), G.probe_vec("gAc"), l, E1, p)          # c in (l), cbar out (l+p)
    B = G.V_cb_c_h(G.probe_vec("gBb"), G.probe_vec("gBc"), l + p, E2, -p)     # c in (l+p), cbar out (l)
    num = G.ghost_line(G.ghost_line(sp.expand(A * B), "gAb", "gBc", 1 / GHOST_NUM), "gBb", "gAc", 1 / GHOST_NUM)
    mine = res4(U.loop_pole(-num, [U.V(), p], 4, 1), "Q05 D03 validation")   # -1 ghost loop
    ref = json.loads((PV / "P16_batchC.json").read_text())["graphs"]["D03"]["result"]
    return {"mine_d4": str(mine), **ratio(mine, hh_reference(ref))}


def D01():
    p, l = U.vec("p1"), U.vec("l")
    E1, E2 = ext_mat("a1"), ext_mat("a2")
    A = G.V3_grav((E1, p), (G.probe_mat("hA1"), -l), (G.probe_mat("hA2"), l - p + 0 * l))
    # lines: hA1 -- hB1 carries l (from A to B: incoming at B = l, at A = -l); hA2 -- hB2 carries (l - p)... at A incoming p + (-l) + x = 0 -> x = l - p
    B = G.V3_grav((E2, -p), (G.probe_mat("hB1"), l), (G.probe_mat("hB2"), p - l))
    num = G.contract_vertices(A, B, [("h", "hA1", "hB1"), ("h", "hA2", "hB2")])
    mine = res4(U.loop_pole(num, [U.V(), -p], 4, sp.Rational(1, 2)), "Q05 D01 validation")
    ref = json.loads((PV / "P16_batchC.json").read_text())["graphs"]["D01"]["result"]
    return {"mine_d4": str(mine), **ratio(mine, hh_reference(ref))}


def D06():
    p, l = U.vec("p1"), U.vec("l")
    gc, gb = U.vec("gc"), U.vec("gb")
    A = G.V_cb_c_h(G.probe_vec("gA"), gc, p, G.probe_mat("hA"), l)            # ext c (p), h in (l), cbar out (p+l)
    B = G.V_cb_c_h(gb, G.probe_vec("gB"), p + l, G.probe_mat("hB"), -l)       # c in (p+l), h in (-l), ext cbar
    num = G.graviton_line(sp.expand(A * B), "hA", "hB", 1)
    num = G.ghost_line(num, "gA", "gB", 1 / GHOST_NUM)
    mine = res4(U.loop_pole(num, [U.V(), p], 4, 1), "Q05 D06 validation")
    D = json.loads((Q / "certificates" / "Q02_bulk_graphs.json").read_text())["D06"]["laurent"]
    H, T = U.H, U.T
    pp = H("p1", "p1")
    ref = sp.sympify(D["S1_delta^rho_nu_p4"]["UV_pole"]) * (H("gb", "gc") + T("gb", "gc")) * pp ** 2 + \
        sp.sympify(D["S2_p^rho_p_nu_p2"]["UV_pole"]) * H("gb", "p1") * H("gc", "p1") * pp + \
        sp.sympify(D["SE_Gt^rho_nu_p4"]["UV_pole"]) * T("gb", "gc") * pp ** 2
    return {"mine_d4": str(mine), **ratio(mine, sp.expand(ref))}


def D13():
    p, l = U.vec("p1"), U.vec("l")
    KE, gc = ext_mat("aK"), U.vec("gc")
    A = G.V_cb_c_h(G.probe_vec("gA"), gc, p, G.probe_mat("hA"), l)            # ext c (p), cbar out (p+l)
    B = G.V_K_c_h(KE, G.probe_vec("gB"), p + l, G.probe_mat("hB"), -l)        # ext K, c in (p+l)
    num = G.graviton_line(sp.expand(A * B), "hA", "hB", 1)
    num = G.ghost_line(num, "gA", "gB", 1 / GHOST_NUM)
    mine = res4(U.loop_pole(num, [U.V(), p], 3, 1), "Q05 D13 validation")
    r = json.loads((PV / "P15_D13.json").read_text())["result"]
    H, T = U.H, U.T
    pp, aKp, gcp = H("p1", "p1"), H("aK", "p1"), H("gc", "p1")
    ref = sp.sympify(r["T2"]["UV_pole"]) * 2 * (H("aK", "gc") + T("aK", "gc")) * aKp * pp + \
        sp.sympify(r["T3"]["UV_pole"]) * aKp ** 2 * gcp + sp.sympify(r["E2"]["UV_pole"]) * 2 * T("aK", "gc") * aKp * pp + \
        sp.sympify(r["T1"]["UV_pole"]) * (H("aK", "aK") + T("aK", "aK")) * gcp * pp
    return {"mine_d4": str(mine), **ratio(mine, sp.expand(ref))}


if __name__ == "__main__":
    name = sys.argv[1]
    t0 = time.time()
    out = {"name": name, **globals()[name](), "seconds": round(time.time() - t0, 1)}
    (Q / "build" / "logs" / f"q05_val_{name}.json").write_text(json.dumps(out, indent=1) + "\n")
    print(name, out["ratio"], out["constant"], out["seconds"])
