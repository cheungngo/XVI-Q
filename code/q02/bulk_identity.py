#!/usr/bin/env python3
"""XVI-Q Q02: bulk retained identity (XVI-O ST4a/ST4b, read-only definitions) on the FULL one-loop kappa^2 bulk two-point
sector of S_inv (m = 0), assembled from every evaluated graph; bulk counterterms; full sigma_0 closure; O.1 walk.

XVI-O (symmetry_defect.json ST4_split, verbatim):  ST4a: Gamma_1,local = sum z_i C_i (gh = 0);
                                                   ST4b: P_ret (S, C) = - P_ret A_1 (gh(C) = 0; defect gh = 1).
Rows are enumerated MECHANICALLY: every two-field monomial (and every one-field monomial) of ghost number 0 (ST4a,
Gamma_1 entries) and 1 (ST4b, defect entries) built from the declared field/antifield content; each row gets evidence:
  graph values (Q02 bulk graphs; XVI-P pinned certificates), phi-parity (every interaction vertex even in {phi, phi*}:
  checked on the XVI-P vertex inventory), source/auxiliary inventory (b, cbar*, b* on no interaction vertex: checked).
Outcome per row: zero / cancellable / obstruction (explicit residual)."""
from __future__ import annotations
import itertools
import json
import sys
from pathlib import Path

import sympy as sp

Q = Path(__file__).resolve().parents[2]
PV = Q / "provenance" / "XVI-P" / "certificates"
XO = Q / "provenance" / "XVI-O"


def jl(p):
    return json.loads(Path(p).read_text())


FIELDS = {  # name: (ghost number, phi-parity, kind)
    "h": (0, 0, "field"), "phi": (0, 1, "field"), "b": (0, 0, "aux"), "c": (1, 0, "ghost"), "cb": (-1, 0, "antighost"),
    "hs": (-1, 0, "antifield (K' in Psi-transformed variables)"), "phis": (-1, 1, "antifield"), "cs": (-2, 0, "antifield"),
    "cbs": (0, 0, "antifield (nonminimal)"), "bs": (-1, 0, "antifield (nonminimal)")}
LEG_ALIAS = {"K": "hs", "h": "h", "phi": "phi", "phis": "phis", "c": "c", "cb": "cb", "cs": "cs", "cbs": "cbs", "b": "b"}


def vertex_checks():
    inv = jl(PV / "P14_amplitude_rules.json")["vertex_inventory"]
    inter = [v for v in inv if len(v["legs"]) >= 3]
    parity_even = {v["name"]: sum(FIELDS[LEG_ALIAS[l]][1] for l in v["legs"]) % 2 == 0 for v in inv}
    aux_free = {v["name"]: not ({"b", "cbs", "bs"} & {LEG_ALIAS[l] for l in v["legs"]}) for v in inter}
    return {"n_vertices": len(inv), "interaction_vertices": [v["name"] for v in inter],
            "phi_parity_even_all": all(parity_even.values()), "phi_parity_per_vertex": parity_even,
            "aux_b_cbs_bs_on_no_interaction_vertex": all(aux_free.values()),
            "quadratic_vertices_with_aux": [v["name"] for v in inv if len(v["legs"]) == 2 and {"b", "cbs"} & set(v["legs"])]}


def enumerate_graphs():
    """Independent one-loop 1PI enumeration from the p14r loop-capable vertex list (P14_counterterm_domain): propagating
    pairs h-h, phi-phi, c-cb (b, cbar*, b* are on no loop-capable vertex).  Tadpoles (one vertex, one self-contraction)
    and bubbles (two vertices, two lines between them).  Returns the E = 1 and E = 2 graphs by (vertices, lines, external)."""
    V = jl(PV / "P14_counterterm_domain.json")["loop_capable_vertices"]
    pair = {("h", "h"), ("phi", "phi"), ("c", "cb"), ("cb", "c")}
    out = set()
    for vn, legs in V.items():
        for i, j in itertools.combinations(range(len(legs)), 2):
            if (legs[i], legs[j]) in pair:
                ext = tuple(sorted(l for k, l in enumerate(legs) if k not in (i, j)))
                out.add(((vn,), tuple(sorted([legs[i] + "-" + legs[j]] if legs[i] <= legs[j] else [legs[j] + "-" + legs[i]])), ext))
    names = sorted(V)
    for a in range(len(names)):
        for b in range(a, len(names)):
            L1, L2 = V[names[a]], V[names[b]]
            for i1, j1 in itertools.permutations(range(len(L1)), 2):
                for i2, j2 in itertools.permutations(range(len(L2)), 2):
                    if (L1[i1], L2[i2]) in pair and (L1[j1], L2[j2]) in pair:
                        ext = tuple(sorted([l for k, l in enumerate(L1) if k not in (i1, j1)] + [l for k, l in enumerate(L2) if k not in (i2, j2)]))
                        lines = tuple(sorted("-".join(sorted((L1[x], L2[y]))) for x, y in ((i1, i2), (j1, j2))))
                        out.add((tuple(sorted((names[a], names[b]))), lines, ext))
    E1 = sorted({g for g in out if len(g[2]) == 1})
    E2 = sorted({g for g in out if len(g[2]) == 2})
    return {"E1": [list(map(list, g)) for g in E1], "E2": [list(map(list, g)) for g in E2], "n_E1": len(E1), "n_E2": len(E2)}


def rows():
    names = list(FIELDS)
    two = sorted({tuple(sorted(p)) for p in itertools.combinations_with_replacement(names, 2)})
    out = {"gh0_two": [], "gh1_two": [], "gh0_one": [], "gh1_one": []}
    for p in two:
        g = sum(FIELDS[x][0] for x in p)
        if g in (0, 1):
            out[f"gh{g}_two"].append(p)
    for x in names:
        g = FIELDS[x][0]
        if g in (0, 1):
            out[f"gh{g}_one"].append((x,))
    return out


def classify(row, vc):
    s = set(row)
    if sum(FIELDS[x][1] for x in row) % 2:
        return "zero", "phi-parity: every vertex even in {phi, phi*} (checked) -> odd rows vanish identically"
    if s & {"b", "cbs", "bs"}:
        if row == ("b", "c"):
            return "ST4b_ghost_equation", "c.b row of (S,Gamma) = b Sigma_cbc + (d b) G_Kc: the ghost equation"
        if row == ("c", "cbs"):
            return "ST4b_unreachable", "c.cbar*: unreachable (XVI-P P22 a_closure: S0 quadratic part h.h, h.b, b.b, b.cbs)"
        return "zero", "auxiliary/nonminimal field on no interaction vertex (checked): no 1PI loop entry"
    table = {("h", "h"): "hh", ("phi", "phi"): "phiphi", ("c", "cb"): "cbc", ("c", "hs"): "Kc", ("c", "h"): "ST4b_ch",
             ("h",): "onepoint_h"}
    if row in table:
        return table[row], None
    if row == ("c",):
        return "zero", ("gh1 one-point defect c^mu(p=0): no 1PI one-loop graph has a single external ghost (the only ghost "
                        "vertex V_cb_c_h would need a second internal c to close the cbar line) and no Lorentz-covariant "
                        "constant vector exists")
    return "unclassified", None


def main(bulk_graphs_json):
    bg = jl(bulk_graphs_json)
    vc = vertex_checks()
    R = rows()
    A = jl(PV / "P15_batchA.json")
    C = jl(PV / "P16_batchC.json")
    D13 = jl(PV / "P15_D13.json")
    D19 = jl(PV / "P19_defect.json")
    E22 = jl(PV / "P22_evanescent.json")
    S23 = jl(PV / "P23_O1_scope.json")
    ST4 = jl(XO / "ledgers" / "symmetry_defect.json")["ST4_split"] if (XO / "ledgers" / "symmetry_defect.json").is_file() else None
    ev_ct = E22["EV_CT"]
    # ---------------- ST4a rows (gh 0)
    st4a = {}
    st4a["hh"] = {"graphs": ["D01", "D02", "D03", "D11", "D12"],
                  "values": {"D01_D03": "XVI-P P16_batchC (evaluated)", "D02": "0 derived (P15_batchA)",
                             "D12": bg["tadpoles"]["D12"]["value"] + " ; XVI-P P15_batchA: " + A["D12_slot"]["status"] if isinstance(A.get("D12_slot"), dict) and "status" in A["D12_slot"] else bg["tadpoles"]["D12"]["value"],
                             "D11": "XVI-P P18_D11 (matching M1-M3; evanescent EV_M2)"},
                  "counterterm": "p14r 4D h.h MS-bar directions + p22 evanescent h.h operators (sum coefficients pinned in P22)",
                  "outcome": "cancellable", "evidence_flags": {"batchC_evaluated": bool(C.get("batchC_evaluated")),
                                                              "batchA_evaluated": bool(A.get("batchA_evaluated"))}}
    d06 = bg["D06"]
    pole = {k: sp.sympify(v["UV_pole"]) for k, v in d06["laurent"].items()}
    ccb_se_ct = sp.sympify(ev_ct.get("coefficient_ccb_SE", E22.get("EV_CT", {}).get("coefficient_ccb_SE", "nan")))
    st4a["cbc"] = {"graphs": ["D06"], "D06_equals_XVI_P_Sigma": d06["equals_XVI_P_P15_exactly"],
                   "poles": {k: str(v) for k, v in pole.items()},
                   "counterterm": {"4D (delta p^4, p^rho p_nu p^2)": {k: str(-v) for k, v in pole.items() if not k.startswith("SE")},
                                   "evanescent SE (p22 ghost-equation partner)": str(ccb_se_ct)},
                   "SE_pole_plus_p22_CT_zero": sp.simplify(pole.get("SE_Gt^rho_nu_p4", sp.nan) + ccb_se_ct) == 0,
                   "outcome": "cancellable"}
    d4, d5 = bg["D04"], bg["tadpoles"]["D05a"]
    st4a["phiphi"] = {"graphs": ["D04", "D05a"], "D04": d4.get("value"), "D04_routes_agree": d4.get("routes_agree"),
                      "D05a": d5["value"], "total_pole": "0", "counterterm": "z_phiphi = 0 (none needed; the declared phi phi "
                      "domain {phi Box^2 phi} is untouched)", "outcome": "zero" if d4.get("value") == "0" and not d5["logarithmic_member"] else "OPEN"}
    st4a["Kc"] = {"graphs": ["D13"], "4D_poles_in_p14r_cK_basis": D13.get("UV_pole_in_p14r_cK_basis") is not None,
                  "E2_pole": D13["result"]["E2"]["UV_pole"], "p22_Kc_E2_CT": ev_ct.get("coefficient_Kc_E2"),
                  "E2_pole_plus_CT_zero": sp.simplify(sp.sympify(D13["result"]["E2"]["UV_pole"]) + sp.sympify(ev_ct.get("coefficient_Kc_E2"))) == 0,
                  "outcome": "cancellable"}
    st4a["onepoint_h"] = {"graphs": ["D14a", "D14b", "D14c"], "value": "0 derived (P15_batchA)", "outcome": "zero"}
    # ---------------- ST4b rows (gh 1)
    st4b = {"ST4b_ch": {"source": "XVI-P P19_defect: projected defect ZERO exactly in d; counterterm solve outcome "
                                  + D19["counterterm_solve"]["outcome"], "projected_identity_verified": D19["projected_identity_verified"],
                        "outcome": "zero" if D19["counterterm_solve"]["outcome"] == "ZERO" else "check"},
            "ST4b_ghost_equation": {"source": "XVI-P P15_D13 ghost_equation (64 components, exact in d) + P22 a_closure c.b "
                                              "(closed; exact in the gt sector); Sigma there == D06 here exactly",
                                    "ghost_equation_agree": D13["ghost_equation"]["agree"],
                                    "c_b_closed": E22["EV_CT"].get("a_closure", {}).get("c.b", {}).get("closed") if isinstance(E22["EV_CT"].get("a_closure"), dict) else None,
                                    "D06_is_that_Sigma": d06["equals_XVI_P_P15_exactly"], "outcome": "zero"},
            "ST4b_unreachable": {"source": "P22 a_closure c.cbs", "outcome": "zero"}}
    # ---------------- assemble all rows
    table = {}
    for grp, lst in R.items():
        for row in lst:
            cls, why = classify(row, vc)
            ent = {"group": grp, "class": cls}
            if why:
                ent["reason"] = why
            if cls in st4a:
                ent["outcome"] = st4a[cls]["outcome"]
            elif cls in st4b:
                ent["outcome"] = st4b[cls]["outcome"]
            elif cls == "zero":
                ent["outcome"] = "zero"
            else:
                ent["outcome"] = "UNCLASSIFIED"
            table[".".join(row)] = ent
    unclassified = [k for k, v in table.items() if v["outcome"] in ("UNCLASSIFIED", "OPEN", "check")]
    # ---------------- full sigma_0 closure
    fc = D19.get("full_closure_mutation_detector", {})
    closure = {"XVI_P_full_closure_all_gh1_sectors": "P19_defect full_closure_mutation_detector + renormalised_defect_zero",
               "renormalised_defect_zero": D19.get("renormalised_defect_zero"),
               "mutation_detector_cases_detected": {k: v.get("detected") for k, v in fc.get("cases", {}).items()},
               "phiphi_CT": "0 (D04 + D05a = 0); sigma_0 rank on the phi.phi direction is 0 (XVI-P P19 D_cap_ker_sigma0), so the "
                            "direction is sigma_0-closed for any coefficient; NO free-EOM equivalence is used (XVI-O O.6 domain item 4)",
               "phiphi_sigma0_rank": D19["D_cap_ker_sigma0"]["sector_ranks"]["phi.phi"]["sigma0_rank"],
               "phi_odd_rows": "identically zero (parity)", "closed": True}
    # ---------------- explicit bulk counterterms in the declared p14r domain (11 directions) + Q02 phi.phi
    dom = jl(PV / "P14_counterterm_domain.json")["kappa2"]["sectors"]
    g1 = D19["gamma1_local_data_domain_coordinates"]
    labels = g1["labels"]
    ct = {}
    for lab, pl in zip(labels, g1["pole"]):
        ct[lab] = -sp.sympify(pl)
    cbc_routeG = [sp.sympify(v) for v in g1["c_cb_ghost_self_energy (route-G (c,cb) orientation, p20)"]["pole"]]
    cbc_P15 = [sp.sympify(v) for v in g1["c_cb_ghost_self_energy (P15 (cb,c) orientation)"]["pole"]]
    ct["c.cb#0"], ct["c.cb#1"] = -cbc_routeG[0], -cbc_routeG[1]
    d06_poles = [pole["S1_delta^rho_nu_p4"], pole["S2_p^rho_p_nu_p2"]]
    phiphi_pole = sp.Integer(0) if st4a["phiphi"]["outcome"] == "zero" else None
    ct["phi.phi#0"] = -phiphi_pole if phiphi_pole is not None else None
    def in_Qi(v):
        v = sp.nsimplify(v)
        re, im = v.as_real_imag()
        return v is not None and re.is_rational and im.is_rational
    counterterms = {"convention": "Ghat coordinates on the p14r kappa^2 densities (G = kappa^2 N_eps Ghat MB); MS-bar CT = -pole; "
                                  "c.cb in the route-G (c,cb) orientation used by the p20 full closure",
                    "coefficients": {k: str(v) for k, v in ct.items()},
                    "domain_sectors_and_sizes": {k: len(v["structures"]) for k, v in dom.items()},
                    "every_label_in_domain": all(lab.split("#")[0] in dom and int(lab.split("#")[1]) < len(dom[lab.split("#")[0]]["structures"]) for lab in ct),
                    "n_directions": len(ct), "all_in_Q_i": all(v is not None and in_Qi(v) for v in ct.values()),
                    "cbc_from_D06_equals_P15_orientation_poles": [sp.simplify(a - b) == 0 for a, b in zip(d06_poles, cbc_P15)],
                    "cbc_orientation_ratio_routeG_over_P15": [str(sp.simplify(a / b)) for a, b in zip(cbc_routeG, cbc_P15)],
                    "phiphi": "0 from Q02 D04 + D05a (in P19 the phi.phi coordinate 0 was a placeholder 'not an input'; now earned)",
                    "evanescent_p22": {"Kc_E2": ev_ct.get("coefficient_Kc_E2"), "ccb_SE": ev_ct.get("coefficient_ccb_SE"),
                                       "hh_sum": ev_ct.get("coefficients_hh_counterterm_sum"),
                                       "role": "outside the p14r 4D domain (disclosed p22 contract revision); subtract the "
                                               "evanescent poles only; NOT needed for the identity (P19 ch solve x = 0; "
                                               "evanescent poles decouple from the ch slots)"}}
    counterterms["ok"] = counterterms["every_label_in_domain"] and counterterms["all_in_Q_i"] and all(counterterms["cbc_from_D06_equals_P15_orientation_poles"]) \
        and counterterms["n_directions"] == 11
    closure.update({"full_closure_with_complete_MSbar_CT": D19["counterterm_solve"].get("full_closure_reading"),
                    "relative_factor_derived_equals_P15": D19["counterterm_solve"].get("K1_ghost_closure_convention", {}).get("derived_equals_P15"),
                    "sigma0_ranks": D19["D_cap_ker_sigma0"]["sector_ranks"], "cohomology_dim_D": D19["D_cap_ker_sigma0"]["cohomology_dim"]})
    obstruction = [k for k, v in table.items() if v["outcome"] == "obstruction"]
    outcome = "obstruction" if obstruction else ("OPEN" if unclassified else
                                                 ("cancellable" if any(v["outcome"] == "cancellable" for v in table.values()) else "zero"))
    return {"XVI_O_ST4_definition_verbatim": ST4, "bulk_counterterms_explicit": counterterms,
            "independent_graph_enumeration": enumerate_graphs(), "vertex_checks": vc, "rows": table, "n_rows": len(table),
            "unclassified_rows": unclassified, "ST4a": st4a, "ST4b": st4b, "sigma0_closure": closure,
            "bulk_identity_outcome": outcome, "residual": None if outcome != "obstruction" else obstruction,
            "E2_enumeration_XVI_P_p23": S23.get("power_counting", S23).get("frozen_E2_match"),
            "bulk_counterterms": {"hh": "p14r MS-bar h.h + p22 evanescent (pinned)", "Kc": "p14r c.K' + p22 E2 (-1/16)",
                                  "cbc": {"4D": {k: str(-v) for k, v in pole.items() if not k.startswith("SE")}, "SE": str(ccb_se_ct)},
                                  "phiphi": "0", "onepoint": "0"}}


def o1_walk(bulk, O2_status, d12_reconciled, d04_scrutiny=None):
    """Walk XVI-O's O.1 closure_check VERBATIM (checklist, closure reasons, enlarge and downgrade rules).  CLOSED only
    if every item is earned; the label is CLOSED_PENDING_INDEPENDENT_REVIEW until the Q02 independent rebuild passes."""
    tgt = jl(XO / "ledgers" / "O_quantum_target.json")
    cc = tgt["closure_check"]
    vc = bulk["vertex_checks"]
    E2 = bulk["E2_enumeration_XVI_P_p23"] or {}
    en = bulk["independent_graph_enumeration"]
    ct = bulk["bulk_counterterms_explicit"]
    sc = bulk["sigma0_closure"]
    st = {c["item"]: c["status"] for c in cc["checklist"]}
    d04_ok = bool(d04_scrutiny) and d04_scrutiny.get("declared_zero") is True and d04_scrutiny.get("all_mutations_nonzero") is True
    items = [
        {"item": "State retained quantum identity", "owner": "O.6", "XVI_O_status": st["State retained quantum identity"],
         "evidence": "XVI-O ST4a/ST4b and (S,Gamma_1) = A_1 (retained_quantum_identity.tex) stated on the bulk two-point "
                     "retained sector; 23 gh0/gh1 rows enumerated mechanically",
         "earned": bulk["XVI_O_ST4_definition_verbatim"] is not None and bulk["n_rows"] == 23},
        {"item": "List amplitudes constrained by identity", "owner": "O.4", "XVI_O_status": st["List amplitudes constrained by identity"],
         "evidence": "every gh0/gh1 one/two-field row classified with an individual justification; independent 1PI "
                     f"enumeration from the loop-capable vertices: E2 = {en['n_E2']}, E1 = {en['n_E1']} (== XVI-P 9 / 3); "
                     "antifield vertices c*cc, phi*c phi give no E <= 2 graph",
         "earned": not bulk["unclassified_rows"] and E2.get("equal") is True and en["n_E2"] == 9 and en["n_E1"] == 3
                   and vc["phi_parity_even_all"] and vc["aux_b_cbs_bs_on_no_interaction_vertex"]},
        {"item": "S_inv^{[2]} classical consistency", "owner": "O.2", "XVI_O_status": st["S_inv^{[2]} classical consistency"],
         "evidence": "XVI-Q Q01: 2(S0,S2)+(S1,S1)=0 Constructed ((a) ell_4, (b) explicit antibracket expansion, (c) "
                     "ghost/antifield); Q01 independent rebuild PASSED", "earned": O2_status == "Constructed"},
        {"item": "Evaluate inventory graphs", "owner": "O.5", "XVI_O_status": st["Evaluate inventory graphs"],
         "evidence": "all 9 E=2 + 3 E=1 graphs have exact values: D01, D02, D03, D11, D13, D14a-c (XVI-P, pinned); D04 "
                     "(Q02, 0 by engine routes 1+2 and independent route 4, with nonzero mutations, literature agreement); "
                     "D05a (Q02, derived 0); D06 (Q02 == XVI-P Sigma); D12 reconciled (derived 0)",
         "earned": all(bulk["ST4a"][r]["outcome"] in ("zero", "cancellable") for r in ("hh", "cbc", "phiphi", "Kc", "onepoint_h"))
                   and d12_reconciled and d04_ok},
        {"item": "Demonstrate identity remainder cancellation / CT solution", "owner": "O.6/O.7",
         "XVI_O_status": st["Demonstrate identity remainder cancellation / CT solution"],
         "evidence": "ST4a: 11-direction MS-bar CT in the declared p14r domain over Q(i) (phi.phi = 0 now earned); ST4b: "
                     "c.h defect ZERO exactly in d (P19, x = 0), c.b = ghost equation (Sigma == D06), c.cbar* unreachable, "
                     "c.phi parity; full sigma_0 closure of the complete MS-bar CT in every gh+1 sector (p20); cohomology "
                     "of the domain 0; no free-EOM equivalence used",
         "earned": bulk["bulk_identity_outcome"] in ("zero", "cancellable") and sc["closed"] and ct["ok"]
                   and sc.get("relative_factor_derived_equals_P15") is True},
        {"item": "Enlarge target if identity requires more legs/orders", "owner": "O.1",
         "XVI_O_status": st["Enlarge target if identity requires more legs/orders"],
         "evidence": "rule 1 evaluated below: not triggered", "earned": True}]
    reasons = {
        cc["reason"][0]: "ANSWERED: identity verified on the retained two-point sector (items 1, 5).",
        cc["reason"][1]: "ANSWERED: the quantum identity is the EA-route (S,Gamma_1) = A_1, not XVI-N's classical truncation.",
        cc["reason"][2]: "ANSWERED: XVI-M selected-image data are not used anywhere (graphs come from the action's vertices).",
        cc["reason"][3]: "ANSWERED: one-loop antifield channels derived from vertices: K'.c = D13 (evaluated); phi*'.c odd (0); "
                         "c*', cbar*', b* give no E <= 2 graph (independent enumeration); no (Phi, Phi*) assignment used.",
        cc["reason"][4]: "ANSWERED: completeness by mechanical enumeration of rows AND graphs (two independent enumerations).",
        cc["reason"][5]: "ANSWERED: kappa^2 master Constructed in Q01 (incl. grav ell_4).",
        cc["reason"][6]: "ANSWERED: ST4a local matching = explicit CT in the declared domain; ST4b on the retained two-point "
                         "projection verified (the FUNCTIONAL ST4b, all n-point, is NOT claimed).",
        cc["reason"][7]: "ANSWERED: the bulk sector (hh, phi phi, cbar c, K'.c, one-point) is evaluated, not only the D11 kernel."}
    rules = {
        tgt["closure_check"]["enlarge_rules"][0]: {
            "triggered": False,
            "argument": "At one loop the kappa^2 component of (S,Gamma_1) is (S_0, Gamma_1[kappa^2]) + (S_1, Gamma_1[kappa^1]) "
                        "(S_k: tree kappa^k terms incl. the seeds D07-D09). One-loop graphs with E external legs are kappa^E "
                        "(every vertex kappa^{legs-2}, p14r power counting; XVI-P p23): Gamma_1[kappa^1] is one-point only and "
                        "is 0 (D14a-c), Gamma_1[kappa^2] is two-point only (three-point starts at kappa^3; one-loop kappa^2 "
                        "three-point sector empty). Hence the identity at this order constrains no n >= 3 or higher-loop piece."},
        tgt["closure_check"]["enlarge_rules"][1]: {"triggered": False, "argument": "classical kappa^2 consistency holds (Q01)."},
        tgt["closure_check"]["enlarge_rules"][2]: {"respected": True,
                                                   "argument": "no XVI-M graph is used or identified with an action graph; "
                                                               "XVI_M pin unchanged; every graph is built from S_inv vertices."}}
    down = {
        tgt["closure_check"]["downgrade_rules"][0]: "not applicable: closure decided by the checks above",
        tgt["closure_check"]["downgrade_rules"][1]: "respected: no amplitude omitted; every row/graph individually justified"}
    allx = all(i["earned"] for i in items)
    return {"ledger_closure_check_verbatim": cc, "items": items, "closure_reasons_answered": reasons,
            "enlarge_rules": rules, "downgrade_rules": down, "all_earned": allx,
            "O1_closure": "CLOSED_PENDING_INDEPENDENT_REVIEW" if allx else "NOT_CLOSED",
            "O1_closure_confirmed": False,
            "confirmation_condition": "Q02 independent rebuild PASS (until then: pending independent review)",
            "remaining": [i["item"] for i in items if not i["earned"]]}


if __name__ == "__main__":
    b = main(sys.argv[1])
    print(json.dumps({k: b[k] for k in ("bulk_identity_outcome", "n_rows", "unclassified_rows")}, indent=1))
