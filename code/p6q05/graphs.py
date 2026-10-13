#!/usr/bin/env python3
"""XVI-Q Q05 P.6 graph UV poles with own rules (gravity.py, validated vs XVI-P D03/D06/D13 (+D01 for V3)).
External legs: leg 1 momentum p1, leg 2 p2, leg 3 p3 = -p1-p2 (all incoming). Gravitons pol a_i (rank-1 sym), ghost c
pol 'gc_i', cbar pol 'gb_i', K pol 'aK' (rank-1 sym), c* pol 'z'. Pole in units of U (as Q04), MS-bar residue at d->4
with hat/tilde kept (T_* = evanescent).  Grassmann orientation: open ghost lines run from the external c to the
external cbar/K/c*; ghost loops carry -1.  usage: graphs.py GID -> build/logs/q05_GID.json"""
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

GH = 1 / (-sp.I)                     # ghost line numerator factor relative to i/k^2 (validated on D03/D06)
p1, p2, l = U.vec("p1"), U.vec("p2"), U.vec("l")
p3 = -(p1 + p2)


def ext(a):
    return TT.sym(U.vec(a), U.vec(a))


def r4(e, g):
    return sp.expand(dred.substitute_d(e, 4, f"Q05 {g}: MS-bar residue"))


def pole(num, shifts, deg, sym, g):
    return r4(U.loop_pole(num, [U.V()] + shifts[1:] if not shifts else shifts, deg, sym), g)


def G08():
    """ghost triangle, external h1 h2 h3, both orientations, ghost loop -1."""
    tot = 0
    for legs in ((("a1", p1), ("a2", p2), ("a3", p3)), (("a1", p1), ("a3", p3), ("a2", p2))):
        (A, kA), (B, kB), (C, kC) = legs
        vA = G.V_cb_c_h(G.probe_vec("xAb"), G.probe_vec("xAc"), l, ext(A), kA)
        vB = G.V_cb_c_h(G.probe_vec("xBb"), G.probe_vec("xBc"), l + kA, ext(B), kB)
        vC = G.V_cb_c_h(G.probe_vec("xCb"), G.probe_vec("xCc"), l + kA + kB, ext(C), kC)
        n = G.ghost_line(sp.expand(vA * vB), "xAb", "xBc", 1 / GH)
        n = G.ghost_line(sp.expand(n * vC), "xBb", "xCc", 1 / GH)
        n = G.ghost_line(n, "xCb", "xAc", 1 / GH)
        tot += U.loop_pole(-n, [U.V(), kA, kA + kB], 4, 1)
    return r4(tot, "G08")


def G12():
    """h(p1) at A (two phi), phi(p2) at B, phi(p3) at C; graviton C->B carries l."""
    vA = G.V_hphiphi(ext("a1"), l + p2, p3 - l)
    vB = G.V_hphiphi(G.probe_mat("hB"), p2, -l - p2)
    vC = G.V_hphiphi(G.probe_mat("hC"), p3, l - p3)
    n = G.graviton_line(sp.expand(vB * vC), "hB", "hC", 1)
    return r4(U.loop_pole(sp.expand(n * vA), [U.V(), p2, -p3], 4, 1), "G12")


def G13():
    """V_h_phi_phi (internal h, ext phi, internal phi) + V_hh_phi_phi (ext h, internal h, ext phi, internal phi);
    both placements of the two external phi."""
    tot = 0
    for pa, pb in ((p2, p3), (p3, p2)):
        vA = G.V_hphiphi(G.probe_mat("hA"), pa, -l - pa)
        vB = G.V_hhphiphi(ext("a1"), G.probe_mat("hB"), pb, l + pa)
        n = G.graviton_line(sp.expand(vA * vB), "hA", "hB", 1)
        tot += U.loop_pole(n, [U.V(), pa], 4, 1)
    return r4(tot, "G13")


def G11():
    """ext c(p1) -> V_cb_c_h (internal h, internal cbar) ; V_h_phi_phi (internal h, ext phi(p2), internal phi) ;
    V_phis_c_phi (ext phi*(p3), internal c, internal phi).  Graviton V1->V2 carries l (incoming at V2)."""
    gc = U.vec("gc1")
    v1 = G.V_cb_c_h(G.probe_vec("xb"), gc, p1, G.probe_mat("h1"), -l)            # cbar out: p1 - l
    v2 = G.V_hphiphi(G.probe_mat("h2"), p2, -l - p2)                               # phi out to V3: l + p2
    v3 = G.V_phis_c_phi(G.probe_vec("xc"), l + p2)                                  # c in: p1 - l ; phi in: l + p2
    n = G.graviton_line(sp.expand(v1 * v2), "h1", "h2", 1)
    n = G.ghost_line(sp.expand(n * v3), "xb", "xc", 1 / GH)
    return r4(U.loop_pole(n, [U.V(), p2, -p1], 3, 1), "G11")


def G11_swapped():
    """G11 with the external phi at p3 and phi* at p2 (same graph, legs exchanged)."""
    gc = U.vec("gc1")
    v1 = G.V_cb_c_h(G.probe_vec("xb"), gc, p1, G.probe_mat("h1"), -l)
    v2 = G.V_hphiphi(G.probe_mat("h2"), p3, -l - p3)
    v3 = G.V_phis_c_phi(G.probe_vec("xc"), l + p3)
    n = G.contract_vertices(v1, v2, [("h", "h1", "h2")])
    n = G.ghost_line(sp.expand(n * v3), "xb", "xc", 1 / GH)
    return r4(U.loop_pole(n, [U.V(), p3, -p1], 3, 1), "G11 swapped")


def G09():
    """ext c(p1) at A, ext cbar(p2) at C, ext h(p3) at B (on the ghost line); internal graviton A--C."""
    gc, gb = U.vec("gc1"), U.vec("gb2")
    vA = G.V_cb_c_h(G.probe_vec("xAb"), gc, p1, G.probe_mat("hA"), l)              # cbar out: p1 + l
    vB = G.V_cb_c_h(G.probe_vec("xBb"), G.probe_vec("xBc"), p1 + l, ext("a3"), p3)  # out: p1 + l + p3
    vC = G.V_cb_c_h(gb, G.probe_vec("xCc"), p1 + l + p3, G.probe_mat("hC"), -l)
    n = G.ghost_line(sp.expand(vA * vB), "xAb", "xBc", 1 / GH)
    n = G.ghost_line(sp.expand(n * vC), "xBb", "xCc", 1 / GH)
    n = G.graviton_line(n, "hA", "hC", 1)
    return r4(U.loop_pole(n, [U.V(), p1, p1 + p3], 4, 1), "G09")


def G07():
    """ext c(p1) at A, ext h(p3) at B (ghost line), ext K(p2) at C (V_K_c_h); internal graviton A--C."""
    gc = U.vec("gc1")
    vA = G.V_cb_c_h(G.probe_vec("xAb"), gc, p1, G.probe_mat("hA"), l)
    vB = G.V_cb_c_h(G.probe_vec("xBb"), G.probe_vec("xBc"), p1 + l, ext("a3"), p3)
    vC = G.V_K_c_h(ext("aK"), G.probe_vec("xCc"), p1 + l + p3, G.probe_mat("hC"), -l)
    n = G.ghost_line(sp.expand(vA * vB), "xAb", "xBc", 1 / GH)
    n = G.ghost_line(sp.expand(n * vC), "xBb", "xCc", 1 / GH)
    n = G.graviton_line(n, "hA", "hC", 1)
    return r4(U.loop_pole(n, [U.V(), p1, p1 + p3], 3, 1), "G07")


def G10():
    """ext c(p1) at A, ext c(p2) at B (each V_cb_c_h), internal graviton A--B, both internal ghosts end at
    V_cs_c_c (ext c*(p3)).  Ordering c*(z) c_A c_B with c_A the line from A."""
    vA = G.V_cb_c_h(G.probe_vec("xAb"), U.vec("gc1"), p1, G.probe_mat("hA"), l)     # out: p1 + l
    vB = G.V_cb_c_h(G.probe_vec("xBb"), U.vec("gc2"), p2, G.probe_mat("hB"), -l)    # out: p2 - l
    vC = G.V_cs_c_c(U.vec("z"), G.probe_vec("xCa"), p1 + l, G.probe_vec("xCb"), p2 - l)
    n = G.graviton_line(sp.expand(vA * vB), "hA", "hB", 1)
    n = G.ghost_line(sp.expand(n * vC), "xAb", "xCa", 1 / GH)
    n = G.ghost_line(n, "xBb", "xCb", 1 / GH)
    return r4(U.loop_pole(n, [U.V(), p1, -p2], 3, 1), "G10")


def _v3(a_ext, k_ext, pa, ka, pb, kb):
    return G.V3_grav((ext(a_ext), k_ext), (G.probe_mat(pa), ka), (G.probe_mat(pb), kb))


def G05():
    """V3 (ext h p1) -- V_h_phi_phi (ext phi p2) -- V_h_phi_phi (ext phi p3); phi line C->B carries l."""
    A = _v3("a1", p1, "hAB", l + p2, "hAC", p3 - l)
    B = G.V_hphiphi(G.probe_mat("hB"), p2, l)
    C = G.V_hphiphi(G.probe_mat("hC"), p3, -l)
    n = G.contract_vertices(A, B, [("h", "hAB", "hB")])
    n = G.contract_vertices(n, C, [("h", "hAC", "hC")])
    return r4(U.loop_pole(n, [U.V(), p2, -p3], 4, 1), "G05")


def G06():
    """V3 (ext h p1) -- V_hh_phi_phi (ext phi p2, phi p3); two graviton lines, symmetry 1/2."""
    A = _v3("a1", p1, "hA1", -l, "hA2", l - p1)
    B = G.V_hhphiphi(G.probe_mat("hB1"), G.probe_mat("hB2"), p2, p3)
    n = G.contract_vertices(A, B, [("h", "hA1", "hB1"), ("h", "hA2", "hB2")])
    return r4(U.loop_pole(n, [U.V(), -p1], 4, sp.Rational(1, 2)), "G06")


def _g34(kslot):
    """V3 (ext h p3) ; B = V_cb_c_h (ext c p1) ; C = V_cb_c_h (ext cbar p2) or V_K_c_h (ext K p2)."""
    A = _v3("a3", p3, "hAB", -l, "hAC", l - p3)
    B = G.V_cb_c_h(G.probe_vec("xb"), U.vec("gc1"), p1, G.probe_mat("hB"), l)
    if kslot:
        C = G.V_K_c_h(ext("aK"), G.probe_vec("xc"), p1 + l, G.probe_mat("hC"), p3 - l)
    else:
        C = G.V_cb_c_h(U.vec("gb2"), G.probe_vec("xc"), p1 + l, G.probe_mat("hC"), p3 - l)
    n = G.contract_vertices(B, C, [("c", "xb", "xc", GH)])
    n = G.contract_vertices(A, n, [("h", "hAB", "hB"), ("h", "hAC", "hC")])
    return n


def G04():
    return r4(U.loop_pole(_g34(False), [U.V(), p1, -p3], 4, 1), "G04")


def G03():
    return r4(U.loop_pole(_g34(True), [U.V(), p1, -p3], 3, 1), "G03")


def G01():
    """V3^3 triangle, one orientation, symmetry 1."""
    A = _v3("a1", p1, "hAB", -(l + p1), "hAC", l)
    B = _v3("a2", p2, "hBA", l + p1, "hBC", -(l + p1 + p2))
    C = _v3("a3", p3, "hCB", l + p1 + p2, "hCA", -l)
    n = G.contract_vertices(A, B, [("h", "hAB", "hBA")])
    n = G.contract_vertices(n, C, [("h", "hBC", "hCB"), ("h", "hAC", "hCA")])
    return r4(U.loop_pole(n, [U.V(), p1, p1 + p2], 4, 1), "G01")


if __name__ == "__main__":
    g = sys.argv[1]
    t0 = time.time()
    res = globals()[g]()
    out = {"graph": "P6-" + g, "pole_over_U_d4_residue": str(res), "n_terms": len(sp.Add.make_args(res)),
           "zero": res == 0, "seconds": round(time.time() - t0, 1)}
    (Q / "build" / "logs" / f"q05_{g}.json").write_text(json.dumps(out, indent=1) + "\n")
    print(g, out["n_terms"], out["zero"], out["seconds"])
