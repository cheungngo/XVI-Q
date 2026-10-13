#!/usr/bin/env python3
"""XVI-Q Q05 certificates (P.6 continued; honest partial):
  Q05_rules_validation.json  own gravity/ghost rules vs XVI-P D01/D03/D06/D13 (constant ratios), two-point ghost equation
  Q05_P6_graphs.json         UV poles of P6-G03..G13 (own rules); G01/G02 Missing
  Q05_P6_rows.json           kappa^3 identity rows at pole level: c.phi.phi (closed, + mutations); c.c.cbar / c.c.K
                             (tree calibrated; loop residual nonzero, UNRESOLVED); h.h.c Missing
  Q05_P6_status.json         outcome per stage, errata, disclosures
  Q05_guard_scan.json        hardened static policy scan (Q04 review)"""
from __future__ import annotations
import hashlib, json, sys, time  # noqa: E401
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
Q = HERE.parents[1]
CERT = Q / "certificates"
LOG = Q / "build" / "logs"
sys.path.insert(0, str(HERE))
GRAPHS = ["G03", "G04", "G05", "G06", "G07", "G08", "G09", "G10", "G11", "G12", "G13"]
DEG3 = {"G03", "G07", "G10", "G11"}


def task_val(name):
    import validate as V
    return name, {"D01": V.D01, "D03": V.D03, "D06": V.D06, "D13": V.D13}[name]()


def task_graph(g):
    import sympy as sp
    import graphs as GR
    t0 = time.time()
    r = {"G03": GR.G03, "G04": GR.G04, "G05": GR.G05, "G06": GR.G06, "G07": GR.G07, "G08": GR.G08, "G09": GR.G09,
         "G10": GR.G10, "G11": GR.G11, "G12": GR.G12, "G13": GR.G13}[g]()
    s = str(r)
    LOG.mkdir(parents=True, exist_ok=True)
    (LOG / f"q05_{g}.json").write_text(json.dumps({"graph": "P6-" + g, "pole_over_U_d4_residue": s}) + "\n")
    return g, {"pole_over_U_d4_residue": s, "sha256": hashlib.sha256(s.encode()).hexdigest(),
               "n_terms": len(sp.Add.make_args(r)), "zero": r == 0,
               "n_evanescent_terms": sum(1 for t in sp.Add.make_args(r) if any(x.name.startswith("T_") for x in t.free_symbols)),
               "momentum_degree": 3 if g in DEG3 else 4, "seconds": round(time.time() - t0, 1)}


def ghost_eq2():
    import sympy as sp
    import rows as R
    import uvpole as U
    import gravity as G
    import tensor as TT
    q, gb = U.vec("p1"), U.vec("gb")
    d06 = R.D06_own("gb", "X", q)
    d06 = R.strip_t(G.contract(sp.expand(d06 * U.dot(G.probe_vec("Y"), U.vec("gc"))), [("X", "Y")]))
    d13 = R.D13_own("gc", q)
    Tm = TT.sym(-q, gb)
    Tm.eta = -sp.Rational(1, 2) * U.dot(-q, gb)
    pair = sp.expand(TT.ev((1, [("m", G.probe_mat("E"), "i", "j"), ("m", Tm, "j", "i")])))
    d13T = R.strip_t(R.contract_id(d13 * sp.I, pair, "K", "E", "h"))
    return str(sp.cancel(d06 / d13T))


def task_ghost_eq3():
    import sympy as sp
    import rows as R
    import uvpole as U
    import tensor as TT
    import graphs as GR
    k, gb = U.V({"p2": 1}), U.vec("gb2")
    Tm = TT.sym(k, gb)
    Tm.eta = -sp.Rational(1, 2) * U.dot(k, gb)
    orig = GR.ext
    GR.ext = lambda a: Tm if a == "aK" else orig(a)
    try:
        g3, g7 = R.strip_t(GR.G03()), R.strip_t(GR.G07())
    finally:
        GR.ext = orig
    g4, g9 = R.strip_t(GR.G04()), R.strip_t(GR.G09())
    return {"ratio_G04_over_G03_K_to_T": str(sp.cancel(g4 / g3)), "ratio_G09_over_G07_K_to_T": str(sp.cancel(g9 / g7)),
            "ratio_total": str(sp.cancel((g4 + g9) / (g3 + g7))),
            "meaning": "ghost equation (Gamma depends on cbar only through K' = K + dPsi/dh): the cbar-c-h graphs equal "
                       "-2 x the K-c-h graphs with K -> T(cbar), graph by graph; same -2 as tree level and two-point"}


def rows_cert():
    import itertools
    import sympy as sp
    import rows as R
    sol, a, b = R.calibrate_cphiphi()
    sig = sol[0]
    main = R.row_cphiphi(sig)
    muts = {"flip_sigma": R.row_cphiphi(-sig)["zero"] is False}
    # drop channels / graphs
    t = R.row_cphiphi(sig)
    muts["drop_S1_Kc_channel"] = t["without_S1_Kc_zero"] is False
    for g in ("G05", "G06", "G12", "G13"):
        keep = [x for x in ("G05", "G06", "G12", "G13") if x != g]
        hpp = sum(R.load(x) for x in keep)
        th = R.strip_t(-sp.I * R.polarise(hpp, "a1", "p1", "gc1"))
        rr = sp.expand(th + R.strip_t(sp.sympify(t["S0_phistar"])) + R.strip_t(sp.sympify(t["S1_Kc"])))
        muts[f"drop_{g}"] = rr != 0
    g05 = R.strip_t(-sp.I * R.polarise(R.load("G05"), "a1", "p1", "gc1"))
    cpp = {"row": "c.phi.phi", "calibration_sigma_phistar": str(sig),
           "tree_identity_closed_with_sigma": sp.expand(a + sig * b) == 0,
           "channels": {k: main[k] for k in ("S0_hphiphi", "S0_phistar", "S1_Kc")}, "residual": main["residual"],
           "closed_at_pole_level": main["zero"], "mutations_detected": muts,
           "G05_s0_projection_zero": g05 == 0,
           "G05_note": "G05's pole is nonzero but its s0-projection vanishes identically, so the row cannot detect "
                       "dropping G05 (recorded, not a failure; G05 enters other projections)",
           "evanescent_note": "every T structure dies under the s0 polarisation (T(p1, c) = 0): the row must vanish "
                              "exactly; it does (no evanescent counterterm needed or available in this row)",
           "ledger_erratum": "Q04 ledger listed only Gamma_phi phi = 0 for (S1, Gamma2); the S1 h-phi-phi vertex x "
                             "Gamma_K c (D13) channel also enters and is REQUIRED (dropping it leaves a nonzero residual)"}
    ghost = {}
    for tg, sc in (("ccb", -1), ("ccK", 1)):
        tree = R.ghost_row(tg, loop=False)
        sig_tree = sp.cancel(-tree["h"] / tree["c"])
        ch, res = R.ghost_row(tg, loop=True, sig_c=sc)
        T = [x for x in res.free_symbols if x.name.startswith("T_")]
        a_ = sp.symbols("x0:4")
        e = sp.expand(ch["h"] * a_[0] + ch["K"] * a_[1] + ch["c1"] * a_[2] + ch["c2"] * a_[3])
        syms = [x for x in e.free_symbols if x not in a_]
        fit = sp.solve(sp.Poly(e, *syms).coeffs(), a_, dict=True)
        nontriv = [f for f in fit if any(v != 0 for v in f.values())]
        ghost["c.c.cbar" if tg == "ccb" else "c.c.K"] = {
            "tree_sigma_c": str(sig_tree), "tree_closes_with_constant_sigma": not sp.sympify(sig_tree).free_symbols,
            "loop_residual_n_terms": len(sp.Add.make_args(res)), "loop_residual_zero": res == 0,
            "loop_residual_4D_part_zero": sp.expand(res.xreplace({x: 0 for x in T})) == 0,
            "loop_residual_sha256": hashlib.sha256(str(res).encode()).hexdigest(),
            "channel_terms": {k: len(sp.Add.make_args(v)) for k, v in ch.items() if k != "c"},
            "any_channel_coefficient_combination_closes": bool(nontriv),
            "status": "OPEN (unresolved): the assembled one-loop row does not close and no rescaling of its four "
                      "channels closes it; NOT reported as an obstruction. Diagnosis so far: the three-point ghost "
                      "equation holds graph by graph (G04 = -2 G03[K->T], G09 = -2 G07[K->T]), so the cbar-c-h / K-c-h "
                      "graphs are mutually consistent; the defect lies in G10 (c* c c; one loop-shift bug already found "
                      "and fixed in Q05), the K / c channels, or the row assembly (Grassmann ordering of the two "
                      "identical c legs). Next step: independent check of G10 and of the c-channel."}
    return {"identity": "(S, Gamma_1^pole) restricted to three-field gh1 monomials (Q04 ledger rows)",
            "normalisation": "all amplitudes 'i x coefficient' / 'i Gamma'; antifield-channel signs calibrated on the "
                             "classical (S0, S1) = 0 (tree, kappa^1; O.2-verified master equation) and used unchanged",
            "c.phi.phi": cpp, "ghost_rows": ghost,
            "h.h.c": {"status": "Missing", "why": "needs G01 (V3^3; evaluation did not finish within the release "
                                                  "budget) and G02 (V4 not constructed)"},
            "others": "h.phi.c, phi.b.c, phi.c.cbs, c.c.phis: zero by phi-parity; aux rows: no loop entry (Q04 ledger)"}


def main():
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=6) as ex:
        fg = [ex.submit(task_graph, g) for g in sorted(GRAPHS, key=lambda g: g not in ("G07", "G08", "G04", "G09"))]
        fv = [ex.submit(task_val, n) for n in ("D01", "D03", "D06", "D13")]
        fe = ex.submit(task_ghost_eq3)
        graphs = dict(f.result() for f in fg)
        vals = dict(f.result() for f in fv)
        geq3 = fe.result()
    val = {"validations": vals,
           "expected_ratios": {"D01": "I", "D03": "I", "D06": "I", "D13": "-I"},
           "all_constant": all(v["constant"] for v in vals.values()),
           "ratios_as_expected": all(vals[k]["ratio"] == r for k, r in {"D01": "I", "D03": "I", "D06": "I", "D13": "-I"}.items()),
           "convention_I": "mine/XVI-P = I is the Q04-established convention factor (own pole = i x XVI-P pole; D11 in Q04)",
           "D13_orientation_disclosure": "D13 (K'.c) matches XVI-P with constant ratio -I on all structures (T1, T2, T3, "
                                         "E2): an overall sign relative to the I convention, i.e. the Grassmann "
                                         "ordering/orientation of the K'.c slot differs between the two codes. Own "
                                         "rows use own D13 consistently; the c.phi.phi row (which uses own D13) closes "
                                         "exactly, so the own orientation is the one consistent with own vertices.",
           "ghost_equation_two_point": {"ratio_D06_over_D13_K_to_T": ghost_eq2(), "expected_constant": True,
                                        "T": "T_mn = sym(k, cbar) - (1/2) eta k.cbar (cbar momentum k)"},
           "ghost_equation_three_point": geq3,
           "propagators": {"graviton": "i/k^2 (I - J/2), Q4S trace 4; derived: EH quadratic gauge invariant, EH + XVI-P "
                                       "gf (+D.D) = k^2 (I - J/2) (XVI-P K = -Box P)", "ghost": "-1/k^2 (from i cbar Box c)"}}
    (CERT / "Q05_rules_validation.json").write_text(json.dumps(val, indent=1) + "\n")
    gr = {"graphs": graphs, "missing": {"G01": "V3^3 triangle: implemented (graphs.G01), evaluation did not finish within "
                                               "the release budget -> Missing (never zero)",
                                        "G02": "needs V4_grav (not constructed in Q05) -> Missing"},
          "carried_from_Q04": {"G14": "evaluated (Q04)", "G15": "evaluated (Q04)", "T1-T3": "derived 0 (scaleless)"},
          "convention": "pole residue at d -> 4 in units of U, hat/tilde kept (T_* evanescent); external momenta p1, p2, "
                        "p3 = -p1-p2; graph-level conventions in code/p6q05/graphs.py docstrings"}
    (CERT / "Q05_P6_graphs.json").write_text(json.dumps(gr, indent=1) + "\n")
    rw = rows_cert()
    (CERT / "Q05_P6_rows.json").write_text(json.dumps(rw, indent=1) + "\n")
    ok_rules = val["all_constant"] and val["ratios_as_expected"]
    cpp_ok = rw["c.phi.phi"]["closed_at_pole_level"] and rw["c.phi.phi"]["G05_s0_projection_zero"] and \
        all(v for k, v in rw["c.phi.phi"]["mutations_detected"].items() if k != "drop_G05")
    st = {"landmark": "XVI-Q", "release": "Q05", "P6_outcome": "Partial",
          "stages": {"rules_V3_ghost_validated": ok_rules,
                     "graphs_evaluated_pole": sorted(graphs) + ["G14", "G15"], "graphs_missing": ["G01", "G02"],
                     "graphs_derived_zero": ["T1", "T2", "T3"],
                     "row_c_phi_phi": "closed at pole level" if cpp_ok else "NOT closed",
                     "row_c_c_cbar": "OPEN (unresolved assembly)", "row_c_c_K": "OPEN (unresolved assembly)",
                     "row_h_h_c": "Missing"},
          "counterterms": {"h_phi_phi + phi*_c_phi sector": "MS-bar CT = -(pole) of G05+G06+G12+G13 and G11 (local; "
                           "Q(i) coefficients as certified); its sigma_0 / (S, CT) consistency in row c.phi.phi is "
                           "exactly the closed row above" if cpp_ok else "Missing",
                           "other sectors": "Missing (rows not closed / graphs missing)",
                           "kappa3_domain_reduction": "NOT done (structural counts only, Q04)"},
          "defect": "c.phi.phi: zero at pole level; c.c.cbar, c.c.K: Open; h.h.c: Missing",
          "sigma0_closure_kappa3": "NOT established (only the c.phi.phi row)",
          "finite_parts": "Missing (pole level only)",
          "O1_closure_unchanged": "CLOSED (scoped)", "bulk_symmetry_restored": False,
          "errata": {"Q05-E1": "Q04 identity ledger row c.phi.phi omitted the (S1 h-phi-phi x Gamma_Kc) channel (medium; "
                               "fixed in Q05 rows)",
                     "Q05-E2": "code/q04/scalar_loop3.py docstring says G13/G14; the graphs are P6-G14/G15 (info; file is "
                               "pinned by the Q04 layer gate and left unchanged)",
                     "Q05-E3": "own D13 orientation = -I x XVI-P (convention, disclosed in Q05_rules_validation)",
                     "Q05-E4": "G10 loop momentum shifts had the wrong sign in the first Q05 evaluation (found while "
                               "assembling the ghost rows; fixed before release; certified G10 is the corrected one)"},
          "review_q04_folded": {"guard_scan": "hardened (Q05_guard_scan: recursive, 7 reviewer forms flagged) AND claim "
                                              "downgraded to best-effort policy scan (not a security boundary)",
                                "gate_unpinned_files": "Q05 gate rejects any code/**/*.py, tools/*.py or tests/*.py not "
                                                       "pinned", "docstring": "Q05-E2"},
          "flags": {"execution_pass": ok_rules, "claim_status": "partial", "scope": "pole level, one loop, kappa^3"}}
    (CERT / "Q05_P6_status.json").write_text(json.dumps(st, indent=1) + "\n")
    import guard_scan2
    (CERT / "Q05_guard_scan.json").write_text(json.dumps(guard_scan2.report(), indent=1) + "\n")
    print(json.dumps(st["stages"]), round(time.time() - t0, 1))


if __name__ == "__main__":
    main()
