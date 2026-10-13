#!/usr/bin/env python3
"""XVI-Q Q06 certificates (P.6 continued; honest partial):
  Q06_ghost_rows.json        c.c.cbar / c.c.K rows: diagnosis + fix (K-channel Grassmann sign), 4D closure, unique fit,
                             mutations, explicit evanescent remainder
  Q06_G10_crosscheck.json    second route for G10 (re-routed loop momentum, swapped Wick order) == Q05 G10
  Q06_V4_validation.json     V4_grav: EH gauge identity n=2 (V3/V2 control) and n=3 (V4/V3) exactly zero; permutation symmetry
  Q06_fastpath.json          ring fast path == Expr path (G05, G06 bit-identical)
  Q06_P6_hhh.json            G01 (V3^3 triangle) and G02 (V3-V4 bubble, 3 channels) UV poles; Bose checks
  Q06_P6_status.json         outcome per stage, errata, open items"""
from __future__ import annotations
import hashlib, json, sys, time  # noqa: E401
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
Q = HERE.parents[1]
CERT = Q / "certificates"
LOG = Q / "build" / "logs"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "p6q05"))


def sha(s):
    return hashlib.sha256(s.encode()).hexdigest()


def task(name):
    import sympy as sp
    t0 = time.time()
    if name in ("G01", "G02_1", "G02_2", "G02_3"):
        import hhh_fast
        o = hhh_fast.run(name)
        o.pop("seconds", None)
        return name, o
    if name == "V4":
        import v4
        r = {"gauge_n2_V3_vs_V2_zero": v4.gauge_check(2) == 0, "gauge_n3_V4_vs_V3_zero": v4.gauge_check(3) == 0}
        r["V4_permutation_symmetric"], r["V4_terms_generic"] = v4.symmetry_check()
        return name, r
    if name == "G10":
        import g10_route2
        import rows as R
        r2 = g10_route2.G10_route2()
        old = R.load("G10")
        return name, {"route2_equals_Q05_G10": sp.expand(r2 - old) == 0, "n_terms": len(sp.Add.make_args(r2)),
                      "sha256_route2": sha(str(r2)), "route": "graviton momentum l - p1; ghost lines l, p1+p2-l; V_c*cc slots swapped (x -1)"}
    if name == "fast":
        import fastpole as FP
        import graphs as GR
        import gravity as G
        import uvpole as U
        import rows as R
        p1, p2, l, p3 = GR.p1, GR.p2, GR.l, GR.p3
        out = {}
        F = FP.Fast(["a1", "a2", "a3", "p1", "p2", "l"] + [f"{x}{y}" for x in ("hA1", "hA2", "hB1", "hB2") for y in "UW"])
        A = F.to(GR._v3("a1", p1, "hA1", -l, "hA2", l - p1))
        B = F.to(G.V_hhphiphi(G.probe_mat("hB1"), G.probe_mat("hB2"), p2, p3))
        r = GR.r4(F.loop_pole(F.contract(A, B, [("h", "hA1", "hB1"), ("h", "hA2", "hB2")], 4), [U.V(), -p1], 4,
                              sp.Rational(1, 2)), "G06 fast")
        out["G06_fast_equals_expr"] = sp.expand(r - R.load("G06")) == 0
        F = FP.Fast(["a1", "a2", "a3", "p1", "p2", "l"] + [f"{x}{y}" for x in ("hAB", "hAC", "hB", "hC") for y in "UW"])
        A = F.to(GR._v3("a1", p1, "hAB", l + p2, "hAC", p3 - l))
        B = F.to(G.V_hphiphi(G.probe_mat("hB"), p2, l))
        C = F.to(G.V_hphiphi(G.probe_mat("hC"), p3, -l))
        n = F.contract(F.contract(A, B, [("h", "hAB", "hB")], 4), C, [("h", "hAC", "hC")], 4)
        r = GR.r4(F.loop_pole(n, [U.V(), p2, -p3], 4, 1), "G05 fast")
        out["G05_fast_equals_expr"] = sp.expand(r - R.load("G05")) == 0
        return name, out
    if name in ("ccK", "ccb"):
        import rows6
        return name, rows6.row(name)
    raise KeyError(name)


def swap(e, mp):
    import sympy as sp
    rep = {}
    for s in e.free_symbols:
        k, x, y = s.name.split("_")
        rep[s] = sp.Symbol("_".join([k] + sorted([mp.get(x, x), mp.get(y, y)])))
    return sp.expand(e.xreplace(rep))


def main():
    import sympy as sp
    CERT.mkdir(exist_ok=True)
    t0 = time.time()
    names = ["G01", "G02_1", "G02_2", "G02_3", "V4", "G10", "fast", "ccK", "ccb"]
    with ProcessPoolExecutor(5) as ex:
        res = dict(ex.map(task, names))
    rows = {"identity": "(S, Gamma_1^pole) restricted to c.c.X monomials, X = cbar / K (Q05 channels)",
            "diagnosis": ("Q05 assembled the K channel (delta S1/delta h|_{X c_j} . Gamma_{K c_i}) with the orientation sign "
                          "of the h channel.  In the K channel c_i sits in Gamma_{K c_i}, to the right of c_j, so reordering to "
                          "c1 c2 costs one extra Grassmann transposition: K_o enters with -sgn(o).  G10 (checked by a second "
                          "route), the c channels and the identical-c antisymmetrisation were correct."),
            "c.c.K": res["ccK"], "c.c.cbar": res["ccb"],
            "status": ("CLOSED on the 4D (non-evanescent) part at pole level for both rows; the remainder is purely evanescent "
                       "(every term carries one T between external polarisations).  Since the row is linear and the 4D part vanishes, local "
                       "evanescent counterterms equal to minus the evanescent pole parts (p22 prescription) remove it; whether those "
                       "counterterms are s0-consistent on the full kappa^3 evanescent domain is OPEN.")}
    (CERT / "Q06_ghost_rows.json").write_text(json.dumps(rows, indent=1) + "\n")
    (CERT / "Q06_G10_crosscheck.json").write_text(json.dumps(res["G10"], indent=1) + "\n")
    (CERT / "Q06_V4_validation.json").write_text(json.dumps(
        {"vertex": "V4_grav = i x delta^4 S_EH (kappa^2 stripped), code/p6q06/v4.py", **res["V4"],
         "XVI-P_comparison": "not possible: no non-scaleless XVI-P graph contains V4 (the V4 tadpole is scaleless); "
                             "validated instead by the exact EH gauge identity delta0 S4 + delta1 S3 = 0"}, indent=1) + "\n")
    (CERT / "Q06_fastpath.json").write_text(json.dumps(res["fast"], indent=1) + "\n")
    g02 = sp.expand(sum(sp.sympify(res[f"G02_{c}"]["pole_over_U_d4_residue"]) for c in "123"))
    g01 = sp.sympify(res["G01"]["pole_over_U_d4_residue"])
    sw = {"a1": "a2", "a2": "a1", "p1": "p2", "p2": "p1"}
    import rows as R

    def b13(e):
        return sp.expand(R.remap(e, {"p1": R.P3, "p2": R.P2, "a1": "a3", "a3": "a1"}) - e) == 0
    hhh = {"graphs": {
        "G01": {"pole_over_U_d4_residue": str(g01), "sha256": sha(str(g01)), "n_terms": len(sp.Add.make_args(g01)),
                "zero": g01 == 0, "bose_12_symmetric": swap(g01, sw) == g01, "bose_13_symmetric": b13(g01), "method": "fastpole ring path, t-truncated"},
        "G02": {"pole_over_U_d4_residue": str(g02), "sha256": sha(str(g02)), "n_terms": len(sp.Add.make_args(g02)),
                "zero": g02 == 0, "bose_12_symmetric": swap(g02, sw) == g02, "bose_13_symmetric": b13(g02),
                "channels_bose_related": swap(sp.sympify(res["G02_1"]["pole_over_U_d4_residue"]), sw)
                == sp.sympify(res["G02_2"]["pole_over_U_d4_residue"]),
                "channel_n_terms": {c: res[f"G02_{c}"]["n_terms"] for c in "123"}, "method": "fastpole ring path, V4 own"}},
        "note": "momentum p3 = -p1-p2 eliminated; Bose symmetry tested under (a1,p1)<->(a2,p2) and (a1,p1)<->(a3,p3) (full S3)"}
    (CERT / "Q06_P6_hhh.json").write_text(json.dumps(hhh, indent=1) + "\n")
    ok_rows = all(res[t]["residual_4D_zero"] and all(res[t]["mutations_detected"].values()) for t in ("ccK", "ccb"))
    st = {"landmark": "XVI-Q", "release": "Q06", "P6_outcome": "Partial",
          "stages": {"ghost_rows_4D_closed": ok_rows, "ghost_rows_evanescent_remainder": "explicit; removable by p22-type evanescent CTs; their kappa^3 s0-consistency OPEN",
                     "G10_second_route": res["G10"]["route2_equals_Q05_G10"],
                     "V4_validated_gauge_identity": res["V4"]["gauge_n3_V4_vs_V3_zero"],
                     "graphs_evaluated_pole": ["G01", "G02", "G03", "G04", "G05", "G06", "G07", "G08", "G09", "G10",
                                               "G11", "G12", "G13", "G14", "G15"],
                     "graphs_missing": [], "graphs_derived_zero": ["T1", "T2", "T3"],
                     "row_c_phi_phi": "closed at pole level (Q05)",
                     "row_c_c_cbar": "closed on 4D part; evanescent remainder OPEN",
                     "row_c_c_K": "closed on 4D part; evanescent remainder OPEN",
                     "row_h_h_c": "Missing (channel ledger recorded; not assembled)"},
          "h_h_c_ledger": ["Gamma_hhh (G01+G02+G08+G14+G15) [h -> s0 h(c)]",
                           "delta S_2/delta h (EH + gauge fixing) . Gamma_{K h c} (G03+G07)",
                           "s1 h = L_c h  . Gamma_hh pole (all kappa^2 hh self-energy graphs)",
                           "delta S_EH3/delta h|_hh . Gamma_{K c} (D13)",
                           "b eliminated: s cbar ~ F(h) => F(h) . Gamma_{cbar c h} (G04+G09)  [channel absent from the Q04 ledger]"],
          "kappa3_CT_domain_reduction": "Missing", "kappa3_defect": "Missing", "kappa3_CT_solve": "Missing",
          "sigma0_closure_kappa3": "NOT established", "finite_parts": "Missing (pole level only)",
          "errata": {"Q06-E1": "Q05 ghost rows reported OPEN because of a Grassmann-sign error in the K-channel assembly "
                               "(rows.ghost_row); not a property of the graphs (low; corrected here, Q05 record kept)",
                     "Q06-E2": "Q04 ledger for h.h.c omits the b-elimination channel F(h).Gamma_{cbar c h} (low; recorded)"},
          "disclosed_revisions": {"Q06-R1": "Q06-R1 (disclosed revision of the frozen Q05 layer): code/p6q05/release_q05.py now treats the pinned Q05_guard_scan n_files_scanned as a LOWER bound (the recursive scan count grows when Q06 adds code files) and excludes Q06-owned files (code/p6q06/, tests/test_q06.py, tools/*q06*) from its unpinned-file check; those are pinned by the Q06 gate. Q05_guard_scan.json regenerated (26 -> 33 files, still clean). release_q05.py is re-pinned by q06_expected.json."},
          "O1_closure_unchanged": "CLOSED (scoped)", "bulk_symmetry_restored": False,
          "flags": {"execution_pass": ok_rows and res["V4"]["gauge_n3_V4_vs_V3_zero"], "claim_status": "partial",
                    "scope": "pole level, one loop, kappa^3"}, "seconds": round(time.time() - t0, 1)}
    (CERT / "Q06_P6_status.json").write_text(json.dumps(st, indent=1) + "\n")
    print("Q06 certs", st["flags"], st["seconds"])


if __name__ == "__main__":
    main()
