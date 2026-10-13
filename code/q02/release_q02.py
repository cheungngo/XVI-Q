#!/usr/bin/env python3
"""XVI-Q Q02 release gate (fails closed; on reject no RELEASE_STATUS.json is left).
(1) the expected-values file code/q02/q02_expected.json is itself pinned (EXPECTED_SHA256 below);
(2) EVERY certified value: each Q01_*/Q02_* certificate (run-specific 'seconds' removed) must equal its pinned copy
    in the expected file -- every number, string and flag, path by path;
(3) semantic checks independent of the pins: per-file intake pins (each True, count == imports.json), the O.2(c)
    ghost/antifield cells (exact cell set, all True, kappa^0..kappa^2, read from the build input P13_master_kappa2 too),
    D04/D05a/D06/D12 values, bulk rows, O.1 walk consistency, hard stops, conventions files pinned."""
from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path

Q = Path(__file__).resolve().parents[2]
CERT = Q / "certificates"
EXPECTED_FILE = Q / "code" / "q02" / "q02_expected.json"
EXPECTED_SHA256 = "5b50d16935d3ff9a8ab8a9eb51ddcaf7d3bbb46870ac6c9208a80e7e0e1e36ec"
STATUS = (Q / "RELEASE_STATUS.json", CERT / "RELEASE_STATUS.json")
if "--layer" in sys.argv:   # Q03: run as an earlier-layer gate; RELEASE_STATUS is owned by the Q03 gate
    STATUS = ()
CELLS = {"af0_grav_Noether", "af0_scalar_Noether", "af0_ghost_gf_transformed", "af1_invariant", "af1_transformed",
         "af2_invariant", "af2_transformed"}
HARD_STOPS = {"bulk_symmetry_restored": False, "insertion": "frozen", "XV_B_opened": False, "XIV_J_opened": False,
              "Route_B_discharged": False, "Missing_never_zero": True}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def strip(o):
    if isinstance(o, dict):
        return {k: strip(v) for k, v in o.items() if k != "seconds"}
    if isinstance(o, list):
        return [strip(v) for v in o]
    return o


def diff_paths(a, b, path=""):
    if isinstance(a, dict) and isinstance(b, dict):
        out = []
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b:
                out.append(f"{path}/{k} (missing)")
            else:
                out += diff_paths(a[k], b[k], f"{path}/{k}")
        return out
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return [f"{path} (length)"]
        return sum((diff_paths(x, y, f"{path}[{i}]") for i, (x, y) in enumerate(zip(a, b))), [])
    return [] if a == b else [path]


def jl(n):
    f = CERT / f"{n}.json"
    return json.loads(f.read_text()) if f.is_file() else None


def checks():
    c = {}
    c["expected_file_pinned"] = EXPECTED_FILE.is_file() and sha(EXPECTED_FILE) == EXPECTED_SHA256
    exp = json.loads(EXPECTED_FILE.read_text()) if EXPECTED_FILE.is_file() else {}
    bad = {}
    for name, ref in exp.get("certificates", {}).items():
        got = jl(name)
        d = ["(certificate missing)"] if got is None else diff_paths(strip(got), ref)
        if d:
            bad[name] = d[:20]
    c["every_certified_value_matches_pinned"] = bool(exp.get("certificates")) and not bad
    c["_value_mismatches"] = bad
    c["conventions_files_pinned"] = all((Q / f).is_file() and sha(Q / f) == s for f, s in exp.get("file_sha256", {}).items()) \
        and len(exp.get("file_sha256", {})) >= 4
    it = jl("Q01_intake") or {}
    imp = json.loads((Q / "intake" / "imports.json").read_text())
    c["per_file_pins_each_true"] = isinstance(it.get("pins"), dict) and len(it["pins"]) == len(imp["files"]) \
        and all(v is True for v in it["pins"].values()) and it.get("n_pins") == len(imp["files"])
    o2 = jl("Q01_O2") or {}
    tbl = o2.get("items", {}).get("c_ghost_antifield", {}).get("kappa2_table", {})
    c["O2c_cells_exact_and_true"] = set(tbl) == CELLS and all(v is True for v in tbl.values())
    m = Q / "build" / "xvip_run" / "certificates" / "P13_master_kappa2.json"
    if m.is_file():
        t = json.loads(m.read_text())["table_by_kappa_and_antifield_number"]
        c["O2c_build_input_all_orders"] = set(t) == {"kappa0", "kappa1", "kappa2"} and all(
            set(t[k]) == CELLS and all(v is True for v in t[k].values()) for k in t)
    else:
        c["O2c_build_input_all_orders"] = False
    c["O2_constructed"] = o2.get("O2_kappa2_master_status") == "Constructed"
    bg = jl("Q02_bulk_graphs") or {}
    d4, d6, tp = bg.get("D04", {}), bg.get("D06", {}), bg.get("tadpoles", {})
    c["D04_zero_two_routes"] = d4.get("value") == "0" and d4.get("routes_agree") is True and d4.get("route1_A") == "0" \
        and d4.get("route2_IBP_A") == "0" and d4.get("N_proportional_to_p2_times_scalar_denominator") is True
    c["D06_equals_XVI_P_and_exact"] = d6.get("equals_XVI_P_P15_exactly") is True and d6.get("fit_exact") is True and \
        {k: v["UV_pole"] for k, v in d6.get("laurent", {}).items()} == {"S1_delta^rho_nu_p4": "-3/8", "S2_p^rho_p_nu_p2": "1/4", "SE_Gt^rho_nu_p4": "-1/8"}
    c["D05a_D12_no_log_member"] = all(tp.get(g, {}).get("logarithmic_member") is False for g in ("D05a", "D12"))
    c["mutations_detect"] = bg.get("mutations", {}).get("D04_drop_trace_part_changes_A") is True and \
        bg.get("mutations", {}).get("D06_S_f_half_changes_coefficients") is True
    rec = jl("Q02_D12_reconciliation") or {}
    c["D12_reconciled"] = rec.get("consistent") is True and rec.get("status") == "Evaluated_vanishing_under_declared_prescription"
    bi = jl("Q02_bulk_identity") or {}
    c["bulk_rows_all_classified"] = bi.get("n_rows") == 23 and bi.get("unclassified_rows") == [] and \
        all(v.get("outcome") in ("zero", "cancellable") for v in bi.get("rows", {}).values())
    c["bulk_outcome_and_closure"] = bi.get("bulk_identity_outcome") == "cancellable" and bi.get("sigma0_closure", {}).get("closed") is True \
        and bi.get("sigma0_closure", {}).get("renormalised_defect_zero") is True
    c["ST4a_pole_ct_cancel"] = bi.get("ST4a", {}).get("cbc", {}).get("SE_pole_plus_p22_CT_zero") is True and \
        bi.get("ST4a", {}).get("Kc", {}).get("E2_pole_plus_CT_zero") is True
    c["vertex_checks"] = bi.get("vertex_checks", {}).get("phi_parity_even_all") is True and \
        bi.get("vertex_checks", {}).get("aux_b_cbs_bs_on_no_interaction_vertex") is True
    o1 = jl("Q02_O1") or {}
    items = o1.get("items", [])
    c["O1_walk_consistent"] = len(items) == 6 and o1.get("O1_closure") == (
        "CLOSED_PENDING_INDEPENDENT_REVIEW" if all(i.get("earned") is True for i in items) else "NOT_CLOSED") \
        and o1.get("all_earned") == all(i.get("earned") is True for i in items) and o1.get("O1_closure_confirmed") is False
    c["O1_enlarge_rule1_not_triggered_and_XVIM_unused"] = list(o1.get("enlarge_rules", {}).values())[0].get("triggered") is False \
        and list(o1.get("enlarge_rules", {}).values())[2].get("respected") is True and len(o1.get("closure_reasons_answered", {})) == 8
    sc = jl("Q02_D04_scrutiny") or {}
    cs = sc.get("cases", {})
    c["D04_scrutiny"] = sc.get("declared_zero") is True and sc.get("all_mutations_nonzero") is True and \
        sc.get("T3_zero_by_Ward_identity") is True and sc.get("N_DRED_equals_2p2(l-p)^2") is True and \
        sc.get("xi_dependence_linear_only (Froeb: m=0 xi^2 terms vanish)") is True and \
        cs.get("MUT non-de Donder gauge T4 piece (beta != -2), g4=1", {}).get("value_over_master") == "pp**2*(D - 1)"
    ctx = bi.get("bulk_counterterms_explicit", {})
    c["bulk_CT_explicit_in_domain"] = ctx.get("ok") is True and ctx.get("n_directions") == 11 and \
        ctx.get("coefficients", {}).get("phi.phi#0") == "0" and ctx.get("coefficients", {}).get("c.cb#0") == "-3/8"
    en = bi.get("independent_graph_enumeration", {})
    c["independent_enumeration_9_3"] = en.get("n_E2") == 9 and en.get("n_E1") == 3
    g = jl("Q01_guard") or {}
    c["guard_18_forms_all_raise"] = g.get("n_forms") == 18 and g.get("missed") == [] and "apart_residue" in g.get("forms", [])
    er_file = json.loads((Q / "intake" / "errata_p23.json").read_text())["errata"]
    c["errata_file_equals_certified"] = er_file == it.get("errata")
    er = {e["id"]: e for e in er_file}
    c["errata_E1_E7_severities"] = list(er) == ["E1", "E2", "E3", "E4", "E5", "E6", "E7"] and er["E2"]["severity"] == "low" \
        and "p23_conventions" in er["E4"]["correction"]
    reg = json.loads((Q / "registration" / "local_registration.json").read_text())
    c["registration_corrected"] = "XVI-A..XVI-I" in reg["how_found"] and "XVI-O has no registration flag" in reg["how_found"]
    c["xvin_label_note"] = "r12" in imp["files"]["provenance/XVI-N/RELEASE_STATUS_r10.json"].get("label_note", "")
    q1 = jl("Q01_release") or {}
    c["Q01_layer_accepted"] = q1.get("publication") == "accepted"
    return c, o1, bi


def main():
    try:
        c, o1, bi = checks()
    except Exception as exc:  # noqa: BLE001
        c, o1, bi = {"error": f"{type(exc).__name__}: {exc}"}, {}, {}
    ok = all(v is True for k, v in c.items() if not k.startswith("_"))
    for s in STATUS:
        s.unlink(missing_ok=True)
    closure = o1.get("O1_closure", "NOT_CLOSED") if ok else "NOT_CLOSED"
    rec = {"landmark": "XVI-Q", "release": "Q02", "version_label": "XVI-Q Q02", "date": "2026-10-10",
           "execution_pass": ok, "publication": "accepted" if ok else "REJECTED", "value_checks": c,
           "O1_closure": closure, "O1_closure_confirmed": False,
           "O1_closure_note": "pending independent review: becomes confirmed only after the Q02 independent rebuild passes",
           "O1_closure_scope": o1.get("scope"), "O1_not_implied": o1.get("not_implied"),
           "O1_checklist": [{k: i[k] for k in ("item", "owner", "earned", "evidence")} for i in o1.get("items", [])],
           "bulk_identity_outcome": bi.get("bulk_identity_outcome"), "bulk_symmetry_restored": False, "symmetry_restored": False,
           "flags": {"bulk_O5_graphs_evaluated": "earned" if ok else "not_earned",
                     "D12_reconciled": "earned" if ok else "not_earned",
                     "bulk_identity_verified_two_point_kappa2": "earned" if ok else "not_earned",
                     "bulk_counterterms_solved": "earned" if ok else "not_earned",
                     "sigma0_closure_bulk": "earned" if ok else "not_earned",
                     "O2_kappa2_master_constructed": "earned (Q01)", "d_guard_q02_verified": "earned" if ok else "not_earned",
                     "author_rebuild_passed": "see reproducibility", "independent_rebuild_passed": "pending",
                     "independent_rebuild_passed_q01": True},
           "hard_stops": HARD_STOPS, "XVI_Q_registered": "locally; not published"}
    cr = CERT / "Q02_release.json"
    if ok:
        txt = json.dumps(rec, indent=1) + "\n"
        for s in STATUS:
            s.write_text(txt)
        cr.write_text(txt)
        (CERT / "Q02_release_rejected.json").unlink(missing_ok=True)
        (Q / "ledgers").mkdir(exist_ok=True)
        (Q / "ledgers" / "Q02_milestones.json").write_text(json.dumps(
            {"landmark": "XVI-Q", "release": "Q02", "O1_closure": closure, "O1_checklist": rec["O1_checklist"],
             "flags": rec["flags"], "hard_stops": HARD_STOPS}, indent=1) + "\n")
    else:
        cr.unlink(missing_ok=True)
        (CERT / "Q02_release_rejected.json").write_text(json.dumps(rec, indent=1, default=str) + "\n")
    print(json.dumps({"publication": rec["publication"], "O1_closure": closure,
                      "failed": [k for k, v in c.items() if not k.startswith("_") and v is not True],
                      "mismatches": c.get("_value_mismatches")}, default=str)[:1500])
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
