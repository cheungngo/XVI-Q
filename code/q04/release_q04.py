#!/usr/bin/env python3
"""XVI-Q Q04 release gate (fails closed; on reject no RELEASE_STATUS.json is left).
(1) code/q04/q04_expected.json pinned; EVERY Q01..Q04 certificate equals its pinned copy; pinned code files re-hashed;
(2) provenance re-hashed (all 33 pinned, no unpinned file); bound Q02 independent report re-hashed;
(3) semantic P.6 checks: inventory == XVI-P 15 + 3 tadpoles; kappa = E; V_hhh_phi_phi two routes + mutation; engine
    validated vs XVI-P D11 (ratio exactly I) + mutation; scalar-loop Ward residue 0 + both mutations nonzero; graph table
    counts; honesty: P6_outcome 'Partial', defect / CT / closure 'Missing' unless every graph is evaluated;
(4) Q01..Q03 layer gates accepted; O.1 record unchanged (CLOSED (scoped), confirmed); every gate read in the source set."""
from __future__ import annotations
import hashlib, json, sys  # noqa: E401
from pathlib import Path

Q = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Q / "tools"))
import source_set as SS  # noqa: E402
CERT = Q / "certificates"
EXPECTED_SHA256 = "9b4e38d5c3894455bd9b20801f6dbf2b1bbfa0c0a46b1e28b446601717514a15"
STATUS = (Q / "RELEASE_STATUS.json", CERT / "RELEASE_STATUS.json")
GENERATED = ("certificates/", "ledgers/", "build/")
HARD_STOPS = {"bulk_symmetry_restored": False, "insertion": "frozen", "XV_B_opened": False, "XIV_J_opened": False,
              "Route_B_discharged": False, "Missing_never_zero": True}
READS: set[str] = set()


def rd(rel):
    READS.add(rel)
    return (Q / rel).read_bytes()


def jr(rel):
    return json.loads(rd(rel))


def sha_b(b):
    return hashlib.sha256(b).hexdigest()


def strip(o):
    if isinstance(o, dict):
        return {k: strip(v) for k, v in o.items() if k != "seconds"}
    if isinstance(o, list):
        return [strip(v) for v in o]
    return o


def _live_scan_clean():
    sys.path.insert(0, str(Q / "code" / "q04"))
    import guard_scan
    READS.add("code/q04/guard_scan.py")
    return guard_scan.scan() == []


def checks():
    c = {}
    eb = rd("code/q04/q04_expected.json") if (Q / "code/q04/q04_expected.json").is_file() else b""
    c["expected_file_pinned"] = bool(eb) and sha_b(eb) == EXPECTED_SHA256
    exp = json.loads(eb) if eb else {}
    bad = [n for n, ref in exp.get("certificates", {}).items()
           if not (CERT / f"{n}.json").is_file() or strip(jr(f"certificates/{n}.json")) != ref]
    c["every_certified_value_matches_pinned"] = len(exp.get("certificates", {})) >= 18 and not bad
    c["_value_mismatches"] = bad
    badf = [f for f, s in exp.get("file_sha256", {}).items() if not (Q / f).is_file() or sha_b(rd(f)) != s]
    c["code_rehashed"] = len(exp.get("file_sha256", {})) >= 15 and not badf
    c["_file_mismatches"] = badf
    imp = jr("intake/imports.json")
    badp = [f for f, m in imp["files"].items() if not (Q / f).is_file() or sha_b(rd(f)) != m["sha256"]]
    unp = [p.relative_to(Q).as_posix() for p in (Q / "provenance").rglob("*") if p.is_file()
           and p.relative_to(Q).as_posix() not in imp["files"]]
    c["provenance_rehashed_all_pinned"] = len(imp["files"]) == 33 and not badp and not unp
    c["bound_q02_report_rehashed"] = sha_b(rd("reproducibility/independent_rebuild_report_q02.json")) == \
        "e4b54702bc768899f22a00a289de0c0aeb2bf9317b87a28c4a462f05ae4bf495"
    c["bound_q03_report_rehashed"] = sha_b(rd("reproducibility/independent_rebuild_report_q03.json")) == "3b73b81f854cb55e6f1fbe2a76f80d7ad8c6e7164915afef2d595c7b36656417"
    fs = exp.get("file_sha256", {})
    c["tools_source_set_pinned"] = all(f in fs for f in ("tools/source_set.py", "tools/build.py", "tools/rebuild_from_source.py"))
    gs = jr("certificates/Q04_guard_scan.json")
    c["guard_leave_sympy_class_disclosed_and_scan_clean"] = gs["scan_clean"] is True and gs["runtime_closable"] is False         and gs["n_files_scanned"] >= 18 and _live_scan_clean() and set(gs["probes"]) >= {"python_term_walk", "str_rename_reparse", "dense_list_eval"}
    inv = jr("certificates/Q04_P6_inventory.json")
    c["inventory_15_plus_3_tadpoles"] = inv["equals_XVI_P_P23_15"] is True and inv["n_E3_loops"] == 15 and \
        inv["n_E3_tadpoles"] == 3 and inv["by_external"] == {"K,c,h": 2, "c,c,cs": 1, "c,cb,h": 2, "c,phi,phis": 1,
                                                              "h,h,h": 7, "h,phi,phi": 5}
    c["kappa_equals_E_and_1pt_2pt_counts"] = inv["power_counting"]["verified_on_all_enumerated_graphs"] is True and \
        inv["n_by_E"].get("1") == 3 and inv["n_by_E"].get("2") == 9
    v = inv["V_hhh_phi_phi"]
    c["vertex_hhhphiphi_two_routes_and_mutation"] = v["two_routes_agree"] is True and v["mutation_tr_h_cubed_detected"] is True
    ev = jr("certificates/Q04_P6_eval.json")
    c["engine_validated_vs_D11"] = ev["validation_D11"]["ratio_mine_over_XVI_P"] == "I" and \
        ev["validation_D11"]["constant_ratio"] is True and ev["validation_mutation_trace_coef_1"]["constant_ratio"] is False
    w = ev["matter_ward_three_point"]
    c["matter_ward_residue_zero_mutations_detect"] = w["declared"]["residue_d4_zero"] is True and \
        w["mutation_wrong_Lie_sign"]["residue_d4_zero"] is False and w["mutation_drop_X2_bubbles"]["residue_d4_zero"] is False
    c["scalar_hhh_bose_symmetric"] = ev["scalar_loop_hhh"]["bose_symmetric_a1_a2"] is True
    tb = ev["graph_table"]
    n_ev = sum(1 for x in tb if x["status"].startswith("UV pole"))
    n_miss = sum(1 for x in tb if x["status"].startswith("Missing"))
    n_zero = sum(1 for x in tb if x["status"].startswith("derived 0"))
    c["graph_table_counts"] = len(tb) == 18 and (n_ev, n_zero, n_miss) == (ev["n_evaluated_pole"], ev["n_derived_zero"], ev["n_missing"]) \
        and n_ev + n_zero + n_miss == 18 and n_zero == 3
    st = jr("certificates/Q04_P6_status.json")
    sg = st["stages"]
    honest = n_miss > 0 and st["P6_outcome"] == "Partial" and sg["defect_assembly"].startswith("Missing") and \
        sg["counterterm_solve_kappa3"] == "Missing" and sg["sigma0_closure_kappa3"] == "Missing"
    c["honest_partial_outcome"] = honest or (n_miss == 0 and st["P6_outcome"] != "Partial")
    c["no_missing_set_to_zero"] = all(not x["status"].startswith("derived 0") or x["id"].startswith("P6-T") for x in tb)
    c["bulk_symmetry_restored_false"] = st["bulk_symmetry_restored"] is False
    for n in ("Q01_release", "Q02_release", "Q03_release"):
        c[f"{n}_layer_accepted"] = jr(f"certificates/{n}.json").get("publication") == "accepted"
    q3 = jr("certificates/Q03_release.json")
    c["O1_record_unchanged"] = q3.get("O1_closure") == "CLOSED (scoped)" and q3.get("O1_closure_confirmed") is True
    src = set(SS.source_files(Q))
    out = sorted(f for f in READS if f not in src and not f.startswith(GENERATED))
    c["every_gate_read_in_source_set"] = not out and len(READS) > 40
    c["_reads_outside_source_set"] = out
    return c, st, q3, ev


def main():
    try:
        c, st, q3, ev = checks()
    except Exception as exc:  # noqa: BLE001
        c, st, q3, ev = {"error": f"{type(exc).__name__}: {exc}"}, {}, {}, {}
    ok = all(v is True for k, v in c.items() if not k.startswith("_"))
    for s in STATUS:
        s.unlink(missing_ok=True)
    out = {"landmark": "XVI-Q", "release": "Q04", "version_label": "XVI-Q Q04", "date": "2026-10-11",
           "kind": "P.6 entry (one-loop kappa^3 three-point): honest partial", "execution_pass": ok,
           "publication": "accepted" if ok else "REJECTED", "value_checks": c,
           "P6_outcome": st.get("P6_outcome") if ok else "NOT_EARNED", "P6_stages": st.get("stages"),
           "P6_rows": st.get("rows_status"),
           "P6_graphs": [{k: x[k] for k in ("id", "external", "status")} for x in ev.get("graph_table", [])],
           "O1_closure": q3.get("O1_closure") if ok else "see Q03", "O1_closure_confirmed": bool(ok and q3.get("O1_closure_confirmed")),
           "O1_closure_scope": q3.get("O1_closure_scope"), "O1_not_implied": q3.get("O1_not_implied"),
           "O1_caveats": q3.get("O1_caveats"), "XVI_O_ledger_O1": "NOT_CLOSED (read-only, unchanged)",
           "bulk_symmetry_restored": False, "symmetry_restored": False,
           "guard_leave_sympy_class": "disclosed; runtime-unclosable; enforced by static scan of package code (clean)",
           "flags": {"P6_inventory": "earned" if ok else "not_earned", "P6_identity_ledger": "earned" if ok else "not_earned",
                     "P6_engine_validated": "earned" if ok else "not_earned",
                     "P6_matter_subsector_ward_pole": "earned" if ok else "not_earned",
                     "P6_evaluated": "partial", "P6_counterterms": "Missing", "P6_sigma0_closure": "Missing",
                     "independent_rebuild_passed_q02": True, "independent_rebuild_passed_q03": True,
                     "independent_rebuild_passed_q04": "pending", "author_rebuild_passed": "see reproducibility"},
           "hard_stops": HARD_STOPS, "XVI_Q_registered": "locally; not published"}
    cr = CERT / "Q04_release.json"
    if ok:
        txt = json.dumps(out, indent=1) + "\n"
        for s in STATUS:
            s.write_text(txt)
        cr.write_text(txt)
        (CERT / "Q04_release_rejected.json").unlink(missing_ok=True)
        (Q / "ledgers" / "Q04_milestones.json").write_text(json.dumps(
            {"landmark": "XVI-Q", "release": "Q04", "P6_outcome": out["P6_outcome"], "flags": out["flags"],
             "hard_stops": HARD_STOPS}, indent=1) + "\n")
    else:
        cr.unlink(missing_ok=True)
        (CERT / "Q04_release_rejected.json").write_text(json.dumps(out, indent=1, default=str) + "\n")
    print(json.dumps({"publication": out["publication"], "P6": out["P6_outcome"],
                      "failed": [k for k, v in c.items() if not k.startswith("_") and v is not True],
                      "mm": {k: v for k, v in c.items() if k.startswith("_") and v}}, default=str)[:1500])
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
