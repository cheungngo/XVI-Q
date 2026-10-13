#!/usr/bin/env python3
"""XVI-Q Q10 release gate (fails closed; on reject no RELEASE_STATUS.json is left).
(1) code/p6q06/q07/q08/q09/q10/q10_expected.json pinned; EVERY Q01..Q10 certificate equals its pinned copy (Q05 guard-scan
    file count: lower bound); EVERY code/**/*.py, tools/*.py, tests/*.py (incl. the Q05 and Q06 gates) pinned and
    re-hashed; an unpinned file is a reject;
(2) live hardened static scan clean;
(3) semantic: h.h.c tree calibrated, 4D residual 0, every mutation detected, unique fit, matter sector exact,
    remainder explicit and evanescent; evanescent-CT naive test recorded as NOT consistent; joint local CT solvable,
    verified, Bose, evanescent, CT shas; status consistent; O.1 unchanged; Q06 layer accepted; Q07-R1 disclosed."""
from __future__ import annotations
import hashlib, json, sys  # noqa: E401
from pathlib import Path

Q = Path(__file__).resolve().parents[6]
sys.path.insert(0, str(Q / "tools"))
import source_set_q06_q07_q08_q09_q10 as SS  # noqa: E402
CERT = Q / "certificates"
EXPECTED_SHA256 = "45e812586e419899e4dd800bb6a7cbdeed8c7ef1462d22f5a467e882fa2ae080"
STATUS = (Q / "RELEASE_STATUS.json", CERT / "RELEASE_STATUS.json")
HARD_STOPS = {"bulk_symmetry_restored": False, "insertion": "frozen", "XV_B_opened": False, "XIV_J_opened": False,
              "Route_B_discharged": False, "Missing_never_zero": True}
READS: set[str] = set()
Q07_CERTS = ["Q07_hhc_row", "Q07_evanescent_CT", "Q07_P6_status", "Q07_release"]
Q08_CERTS = ["Q08_ghost_equation", "Q08_free_directions", "Q08_P6_status", "Q08_release"]
Q09_CERTS = ["Q09_guard_scan3", "Q09_q04_report_binding", "Q09_trivial_families", "Q09_P6_status"]
Q10_CERTS = ["Q10_family_anatomy", "Q10_erratum_Q08_counts", "Q10_new_manuscripts", "Q10_P6_status"]
Q04_REPORT_SHA256 = "bc2b16d99f01cf172a5c67f6ee1abb5fdb7b9b5c30e6ad3e510c241ffa442a45"
SELF = "code/p6q06/q07/q08/q09/q10/release_q10.py"


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


def code_files():
    out = {p.relative_to(Q).as_posix() for p in (Q / "code").rglob("*.py") if "__pycache__" not in p.parts}
    out |= {p.relative_to(Q).as_posix() for p in (Q / "tools").glob("*.py")}
    out |= {p.relative_to(Q).as_posix() for p in (Q / "tests").glob("*.py")}
    # Q11-R1 (disclosed scope revision): Q11 files (code/p6q06/q07/q08/q09/q10/q11/, tools/*q11*) are pinned by the Q11 gate
    out = {f for f in out if not (f.startswith("code/p6q06/q07/q08/q09/q10/q11/") or (f.startswith("tools/") and "q11" in f))}
    return out - {SELF}


def checks():
    c = {}
    ep = "code/p6q06/q07/q08/q09/q10/q10_expected.json"
    eb = rd(ep) if (Q / ep).is_file() else b""
    c["expected_file_pinned"] = bool(eb) and sha_b(eb) == EXPECTED_SHA256
    exp = json.loads(eb) if eb else {}
    certs = exp.get("certificates", {})

    def _layer(n, o):
        if n in ("Q05_guard_scan", "Q09_guard_scan3") and isinstance(o, dict):  # file counts are lower bounds (Q11-R1)
            o = dict(o)
            ok_n = isinstance(o.get("n_files_scanned"), int) and o["n_files_scanned"] >= certs[n]["n_files_scanned"]
            o["n_files_scanned"] = certs[n]["n_files_scanned"] if ok_n else -1
        return o
    bad = [n for n, ref in certs.items() if not (CERT / f"{n}.json").is_file()
           or _layer(n, strip(jr(f"certificates/{n}.json"))) != ref]
    c["every_certified_value_matches_pinned"] = len(certs) >= 44 and set(Q07_CERTS + Q08_CERTS + Q09_CERTS + Q10_CERTS + ["Q09_release"]) <= set(certs) and not bad
    fs = exp.get("file_sha256", {})
    c["pinned_files_rehashed"] = bool(fs) and all((Q / f).is_file() and sha_b(rd(f)) == h for f, h in fs.items())
    c["no_unpinned_code_tools_tests_file"] = bool(fs) and code_files() <= set(fs) and \
        {"code/p6q05/release_q05.py", "code/p6q06/release_q06.py", "code/p6q06/q07/release_q07.py", "code/p6q06/q07/q08/release_q08.py", "code/p6q06/q07/q08/q09/release_q09.py"} <= set(fs)
    sys.path.insert(0, str(Q / "code" / "p6q05"))
    import guard_scan2
    READS.add("code/p6q05/guard_scan2.py")
    live = guard_scan2.report()
    c["scan_live_clean_incl_q07"] = live["scan_clean"] is True and live["all_reviewer_evasions_flagged"] is True
    sys.path.insert(0, str(Q / "code/p6q06/q07/q08/q09"))  # scan3 module (pinned)
    import guard_scan3
    READS.add("code/p6q06/q07/q08/q09/guard_scan3.py")
    READS.add("code/p6q06/q07/q08/q09/guard_scan3_rules.json")
    s3 = guard_scan3.report()
    c["scan3_live_clean_probes_flagged"] = s3["extension_clean"] is True and s3["all_probes_flagged"] is True and \
        s3["base_scan2_clean"] is True and len(s3["probes_flagged"]) >= 10
    c["Q04_report_bound_by_hash"] = sha_b(rd("reproducibility/independent_rebuild_report_q04.json")) == Q04_REPORT_SHA256
    a = jr("certificates/Q10_family_anatomy.json")
    fm = a["families"]
    c["anatomy_values"] = a["blind_coordinates"]["n"] == 42 and a["blind_coordinates"]["by_sector"] == \
        {"K": 30, "C": 0, "S": 12, "H": 0} and a["blind_coordinates"]["dim_inside_ghost_eq"] == 24 and \
        fm["F"]["n"] == 204 and fm["F"]["zero_vector"] == 15 and fm["F"]["in_N"] == 124 and \
        fm["F"]["lam_any"] == 124 and fm["F"]["in_N_blind"] == 62 and fm["F"]["lam_none"] == 64 and \
        fm["G"]["lam_none"] == 28 and fm["G"]["lam_unique"] == {"I/2": 5} and fm["P"]["lam_none"] == 132 and \
        fm["P"]["lam_unique"] == {"-2": 2} and a["dim_N"] == 152 and a["dim_N_with_ghost_eq"] == 104
    e = jr("certificates/Q10_erratum_Q08_counts.json")
    c["erratum_withdraws_counts"] = e["id"] == "Q10-E1" and len(e["withdrawn"]["Q08"]) == 3 and \
        len(e["withdrawn"]["Q09"]) == 2 and e["restated_bound"]["Q08_value_recorded"] == 36 and \
        "NOT verified" in e["restated_bound"]["caveat"]
    nm = jr("certificates/Q10_new_manuscripts.json")
    c["new_manuscripts_context_only"] = nm["premises_imported"] == [] and nm["XV_B_opened"] is False and \
        nm["sha256"] == sha_b(rd("intake/new_manuscripts_q10.json"))
    st = jr("certificates/Q10_P6_status.json")
    gl = st["goals"]
    c["status_scoped_honest"] = st["P6_outcome"] == "Partial" and st["flags"]["claim_status"] == "partial" and \
        st["sigma0_status_kappa3"] == ("evanescent kappa^3 three-point defect sigma0-exact at pole level; X non-unique "
                                       "(ghost-eq compatible)") and all(gl[k].startswith("OPEN") for k in (
            "F_trace_failures", "G_P_families", "classification_remaining", "invariant_evanescent_basis",
            "p22_kappa2_consistency")) and "Q10-E1" in st["errata"] and \
        "any count of BRST-trivial directions" in st["not_implied"]
    c["revision_disclosed"] = "Q10-R1" in st["disclosed_revisions"]
    q9r = jr("certificates/Q09_release.json") if (CERT / "Q09_release.json").is_file() else {}
    c["Q09_layer_accepted"] = q9r.get("decision") == "ACCEPT"
    q8 = jr("certificates/Q08_P6_status.json")
    c["Q08_result_unchanged"] = q8["goals"]["1_ghost_equation_compatibility"].startswith("ESTABLISHED")
    o1 = jr("certificates/Q03_O1_record.json")
    c["O1_record_unchanged"] = o1["O1_closure"] == "CLOSED (scoped)" and o1["O1_closure_confirmed"] is True and \
        st["O1_closure_unchanged"] == "CLOSED (scoped)" and st["bulk_symmetry_restored"] is False
    c["gate_reads_in_source_set"] = READS <= set(SS.source_files(Q)) | {f"certificates/{n}.json" for n in certs} | \
        {"certificates/Q03_O1_record.json", "certificates/Q09_release.json", "intake/new_manuscripts_q10.json"}
    return c, st


def main():
    for p in STATUS:
        p.unlink(missing_ok=True)
    try:
        c, st = checks()
    except Exception as e:  # noqa: BLE001
        c, st = {"exception": False, "error": repr(e)}, {}
    ok = all(v is True for k, v in c.items() if k != "error")
    rec = {"landmark": "XVI-Q", "release": "Q10", "decision": "ACCEPT" if ok else "REJECT", "checks": c}
    (CERT / "Q10_release.json").write_text(json.dumps(rec, indent=1) + "\n")
    if not ok:
        print(json.dumps(rec, indent=1))
        return 1
    status = {"landmark": "XVI-Q", "release": "Q10", "status": "released (honest partial: erratum Q10-E1 on trivial counts, family anatomy)",
              "P6_outcome": st["P6_outcome"], "sigma0_status_kappa3": st["sigma0_status_kappa3"],
              "O1_closure": "CLOSED (scoped)", "O1_closure_confirmed": True, "hard_stops": HARD_STOPS,
              "independent_rebuild_q05": "not requested", "independent_rebuild_q06": "not requested",
              "independent_rebuild_q07": "PASSED (reviewer; report sha 03f5393f)", "independent_rebuild_q08": "PASSED on execution (reviewer; report sha cc67001e)", "independent_rebuild_q09": "not requested", "independent_rebuild_q10": "pending", "goals": st["goals"], "execution_pass": True, "claim_status": "partial", "gate": SELF}
    for p in STATUS:
        p.write_text(json.dumps(status, indent=1) + "\n")
    (Q / "ledgers").mkdir(exist_ok=True)
    (Q / "ledgers" / "Q10_milestones.json").write_text(json.dumps({"landmark": "XVI-Q", "release": "Q10", **strip(st)},
                                                                  indent=1) + "\n")
    print("ACCEPT")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
