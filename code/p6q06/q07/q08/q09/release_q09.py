#!/usr/bin/env python3
"""XVI-Q Q09 release gate (fails closed; on reject no RELEASE_STATUS.json is left).
(1) code/p6q06/q07/q08/q09/q09_expected.json pinned; EVERY Q01..Q09 certificate equals its pinned copy (Q05 guard-scan
    file count: lower bound); EVERY code/**/*.py, tools/*.py, tests/*.py (incl. the Q05 and Q06 gates) pinned and
    re-hashed; an unpinned file is a reject;
(2) live hardened static scan clean;
(3) semantic: h.h.c tree calibrated, 4D residual 0, every mutation detected, unique fit, matter sector exact,
    remainder explicit and evanescent; evanescent-CT naive test recorded as NOT consistent; joint local CT solvable,
    verified, Bose, evanescent, CT shas; status consistent; O.1 unchanged; Q06 layer accepted; Q07-R1 disclosed."""
from __future__ import annotations
import hashlib, json, sys  # noqa: E401
from pathlib import Path

Q = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(Q / "tools"))
import source_set_q06_q07_q08_q09 as SS  # noqa: E402
CERT = Q / "certificates"
EXPECTED_SHA256 = "b998ce8dc9a332093c9bd0f04bab8a32bc2cbe31efb7d23a303240733bef110e"
STATUS = (Q / "RELEASE_STATUS.json", CERT / "RELEASE_STATUS.json")
HARD_STOPS = {"bulk_symmetry_restored": False, "insertion": "frozen", "XV_B_opened": False, "XIV_J_opened": False,
              "Route_B_discharged": False, "Missing_never_zero": True}
READS: set[str] = set()
Q07_CERTS = ["Q07_hhc_row", "Q07_evanescent_CT", "Q07_P6_status", "Q07_release"]
Q08_CERTS = ["Q08_ghost_equation", "Q08_free_directions", "Q08_P6_status", "Q08_release"]
Q09_CERTS = ["Q09_guard_scan3", "Q09_q04_report_binding", "Q09_trivial_families", "Q09_P6_status"]
Q04_REPORT_SHA256 = "bc2b16d99f01cf172a5c67f6ee1abb5fdb7b9b5c30e6ad3e510c241ffa442a45"
SELF = "code/p6q06/q07/q08/q09/release_q09.py"


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
    # Q10-R1 (disclosed scope revision): Q10 files (code/p6q06/q07/q08/q09/q10/, tools/*q10*) are pinned by the Q10 gate
    out = {f for f in out if not (f.startswith("code/p6q06/q07/q08/q09/q10/") or (f.startswith("tools/") and "q10" in f))}
    return out - {SELF}


def checks():
    c = {}
    ep = "code/p6q06/q07/q08/q09/q09_expected.json"
    eb = rd(ep) if (Q / ep).is_file() else b""
    c["expected_file_pinned"] = bool(eb) and sha_b(eb) == EXPECTED_SHA256
    exp = json.loads(eb) if eb else {}
    certs = exp.get("certificates", {})

    def _layer(n, o):
        if n in ("Q05_guard_scan", "Q09_guard_scan3") and isinstance(o, dict):  # Q10-R1: file counts are lower bounds
            o = dict(o)
            ok_n = isinstance(o.get("n_files_scanned"), int) and o["n_files_scanned"] >= certs[n]["n_files_scanned"]
            o["n_files_scanned"] = certs[n]["n_files_scanned"] if ok_n else -1
        return o
    bad = [n for n, ref in certs.items() if not (CERT / f"{n}.json").is_file()
           or _layer(n, strip(jr(f"certificates/{n}.json"))) != ref]
    c["every_certified_value_matches_pinned"] = len(certs) >= 44 and set(Q07_CERTS + Q08_CERTS + Q09_CERTS) <= set(certs) and not bad
    fs = exp.get("file_sha256", {})
    c["pinned_files_rehashed"] = bool(fs) and all((Q / f).is_file() and sha_b(rd(f)) == h for f, h in fs.items())
    c["no_unpinned_code_tools_tests_file"] = bool(fs) and code_files() <= set(fs) and \
        {"code/p6q05/release_q05.py", "code/p6q06/release_q06.py", "code/p6q06/q07/release_q07.py", "code/p6q06/q07/q08/release_q08.py"} <= set(fs)
    sys.path.insert(0, str(Q / "code" / "p6q05"))
    import guard_scan2
    READS.add("code/p6q05/guard_scan2.py")
    live = guard_scan2.report()
    c["scan_live_clean_incl_q07"] = live["scan_clean"] is True and live["all_reviewer_evasions_flagged"] is True
    sys.path.insert(0, str(Q / "code/p6q06/q07/q08/q09"))
    import guard_scan3
    READS.add("code/p6q06/q07/q08/q09/guard_scan3.py")
    READS.add("code/p6q06/q07/q08/q09/guard_scan3_rules.json")
    s3 = guard_scan3.report()
    c["scan3_live_clean_probes_flagged"] = s3["extension_clean"] is True and s3["all_probes_flagged"] is True and \
        s3["base_scan2_clean"] is True and len(s3["probes_flagged"]) >= 10
    c["Q04_report_bound_by_hash"] = sha_b(rd("reproducibility/independent_rebuild_report_q04.json")) == Q04_REPORT_SHA256
    t = jr("certificates/Q09_trivial_families.json")
    fm = t["families"]
    c["derived_lambdas_honest"] = fm["F"]["lambda_derived"] == "-1" and fm["P"]["lambda_derived"] == "-1" and \
        fm["G"]["lambda_derived"] == "1" and fm["F"]["n_in_N"] == 139 and fm["F"]["n_fail"] == 65 and \
        fm["F"]["n_fail_with_trace_monomials"] == 65 and fm["G"]["n_in_N"] == 0 and fm["P"]["n_in_N"] == 0 and \
        t["dim_trivial_span_derived"] == 59 and "OPEN" in t["diagnosis"]
    st = jr("certificates/Q09_P6_status.json")
    gl = st["goals"]
    c["status_scoped_honest"] = st["P6_outcome"] == "Partial" and st["flags"]["claim_status"] == "partial" and \
        st["sigma0_status_kappa3"] == ("evanescent kappa^3 three-point defect sigma0-exact at pole level; X non-unique "
                                       "(ghost-eq compatible)") and gl["c_invariant_operator_basis"].startswith("OPEN") \
        and gl["d_p22_consistency"].startswith("OPEN") and gl["a_derived_normalisations"].startswith("Partial") and \
        gl["b_classification"].startswith("Partial") and all(v is True for v in st["review_fixes"].values())
    c["revision_disclosed"] = "Q09-R1" in st["disclosed_revisions"]
    q8 = jr("certificates/Q08_P6_status.json")
    c["Q08_result_unchanged"] = q8["goals"]["1_ghost_equation_compatibility"].startswith("ESTABLISHED")
    o1 = jr("certificates/Q03_O1_record.json")
    c["O1_record_unchanged"] = o1["O1_closure"] == "CLOSED (scoped)" and o1["O1_closure_confirmed"] is True and \
        st["O1_closure_unchanged"] == "CLOSED (scoped)" and st["bulk_symmetry_restored"] is False
    q8r = jr("certificates/Q08_release.json") if (CERT / "Q08_release.json").is_file() else {}
    c["Q08_layer_accepted"] = q8r.get("decision") == "ACCEPT"
    c["gate_reads_in_source_set"] = READS <= set(SS.source_files(Q)) | {f"certificates/{n}.json" for n in certs} | \
        {"certificates/Q03_O1_record.json", "certificates/Q08_release.json"}
    return c, st


def main():
    for p in STATUS:
        p.unlink(missing_ok=True)
    try:
        c, st = checks()
    except Exception as e:  # noqa: BLE001
        c, st = {"exception": False, "error": repr(e)}, {}
    ok = all(v is True for k, v in c.items() if k != "error")
    rec = {"landmark": "XVI-Q", "release": "Q09", "decision": "ACCEPT" if ok else "REJECT", "checks": c}
    (CERT / "Q09_release.json").write_text(json.dumps(rec, indent=1) + "\n")
    if not ok:
        print(json.dumps(rec, indent=1))
        return 1
    status = {"landmark": "XVI-Q", "release": "Q09", "status": "released (honest partial: P.6 evanescent CT structure, review fixes)",
              "P6_outcome": st["P6_outcome"], "sigma0_status_kappa3": st["sigma0_status_kappa3"],
              "O1_closure": "CLOSED (scoped)", "O1_closure_confirmed": True, "hard_stops": HARD_STOPS,
              "independent_rebuild_q05": "not requested", "independent_rebuild_q06": "not requested",
              "independent_rebuild_q07": "PASSED (reviewer; report sha 03f5393f)", "independent_rebuild_q08": "in progress (reviewer)", "independent_rebuild_q09": "pending", "goals": st["goals"], "execution_pass": True, "claim_status": "partial", "gate": SELF}
    for p in STATUS:
        p.write_text(json.dumps(status, indent=1) + "\n")
    (Q / "ledgers").mkdir(exist_ok=True)
    (Q / "ledgers" / "Q09_milestones.json").write_text(json.dumps({"landmark": "XVI-Q", "release": "Q09", **strip(st)},
                                                                  indent=1) + "\n")
    print("ACCEPT")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
