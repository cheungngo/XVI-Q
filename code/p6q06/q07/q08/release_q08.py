#!/usr/bin/env python3
"""XVI-Q Q08 release gate (fails closed; on reject no RELEASE_STATUS.json is left).
(1) code/p6q06/q07/q08/q08_expected.json pinned; EVERY Q01..Q08 certificate equals its pinned copy (Q05 guard-scan
    file count: lower bound); EVERY code/**/*.py, tools/*.py, tests/*.py (incl. the Q05 and Q06 gates) pinned and
    re-hashed; an unpinned file is a reject;
(2) live hardened static scan clean;
(3) semantic: h.h.c tree calibrated, 4D residual 0, every mutation detected, unique fit, matter sector exact,
    remainder explicit and evanescent; evanescent-CT naive test recorded as NOT consistent; joint local CT solvable,
    verified, Bose, evanescent, CT shas; status consistent; O.1 unchanged; Q06 layer accepted; Q07-R1 disclosed."""
from __future__ import annotations
import hashlib, json, sys  # noqa: E401
from pathlib import Path

Q = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(Q / "tools"))
import source_set_q06_q07_q08 as SS  # noqa: E402
CERT = Q / "certificates"
EXPECTED_SHA256 = "7e24cb02449aefaabad1fe2e849a722217d93e12290fd8ee0a243403d74996d5"
STATUS = (Q / "RELEASE_STATUS.json", CERT / "RELEASE_STATUS.json")
HARD_STOPS = {"bulk_symmetry_restored": False, "insertion": "frozen", "XV_B_opened": False, "XIV_J_opened": False,
              "Route_B_discharged": False, "Missing_never_zero": True}
READS: set[str] = set()
Q07_CERTS = ["Q07_hhc_row", "Q07_evanescent_CT", "Q07_P6_status", "Q07_release"]
Q08_CERTS = ["Q08_ghost_equation", "Q08_free_directions", "Q08_P6_status"]
SELF = "code/p6q06/q07/q08/release_q08.py"


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
    # Q09-R1 (disclosed scope revision): Q09 files (code/p6q06/q07/q08/q09/, tools/*q09*) are pinned by the Q09 gate
    out = {f for f in out if not (f.startswith("code/p6q06/q07/q08/q09/") or (f.startswith("tools/") and "q09" in f))}
    return out - {SELF}


def checks():
    c = {}
    ep = "code/p6q06/q07/q08/q08_expected.json"
    eb = rd(ep) if (Q / ep).is_file() else b""
    c["expected_file_pinned"] = bool(eb) and sha_b(eb) == EXPECTED_SHA256
    exp = json.loads(eb) if eb else {}
    certs = exp.get("certificates", {})

    def _layer(n, o):
        if n == "Q05_guard_scan" and isinstance(o, dict):
            o = dict(o)
            ok_n = isinstance(o.get("n_files_scanned"), int) and o["n_files_scanned"] >= certs[n]["n_files_scanned"]
            o["n_files_scanned"] = certs[n]["n_files_scanned"] if ok_n else -1
        return o
    bad = [n for n, ref in certs.items() if not (CERT / f"{n}.json").is_file()
           or _layer(n, strip(jr(f"certificates/{n}.json"))) != ref]
    c["every_certified_value_matches_pinned"] = len(certs) >= 39 and set(Q07_CERTS + Q08_CERTS) <= set(certs) and not bad
    fs = exp.get("file_sha256", {})
    c["pinned_files_rehashed"] = bool(fs) and all((Q / f).is_file() and sha_b(rd(f)) == h for f, h in fs.items())
    c["no_unpinned_code_tools_tests_file"] = bool(fs) and code_files() <= set(fs) and \
        {"code/p6q05/release_q05.py", "code/p6q06/release_q06.py", "code/p6q06/q07/release_q07.py"} <= set(fs)
    sys.path.insert(0, str(Q / "code" / "p6q05"))
    import guard_scan2
    READS.add("code/p6q05/guard_scan2.py")
    live = guard_scan2.report()
    c["scan_live_clean_incl_q07"] = live["scan_clean"] is True and live["all_reviewer_evasions_flagged"] is True
    g = jr("certificates/Q08_ghost_equation.json")
    pt = g["particular_solution"] or {}
    c["ghost_equation_compatible"] = g["loops_satisfy_with_rho_minus_2"] is True and \
        g["loop_poles"] == {"rho_total": "-2", "rho_G04_G03": "-2", "rho_G09_G07": "-2"} and \
        g["solvable_with_ghost_equation"] is True and pt.get("verified_A_x_eq_b") is True and \
        g["dim_free_without_ghost_eq"] == 152 and g["dim_free_with_ghost_eq"] == 104
    f = jr("certificates/Q08_free_directions.json")
    c["free_directions_not_overclaimed"] = f["homogeneous_space"]["dim_N"] == 152 and \
        f["dim_N_cap_all_families_best_lambda"] <= f["dim_N_cap_all_families_upper_bound"] < 152 and \
        "FITTED" in f["reading"] and "not shown trivial" in f["reading"]
    st = jr("certificates/Q08_P6_status.json")
    gl = st["goals"]
    c["status_scoped_honest"] = st["P6_outcome"] == "Partial" and st["flags"]["claim_status"] == "partial" and \
        gl["2_invariant_operator_basis"].startswith("Partial") and gl["3_free_parameters"].startswith("Partial") and \
        gl["4_p22_consistency"].startswith("OPEN") and gl["1_ghost_equation_compatibility"].startswith("ESTABLISHED") \
        and "uniqueness of X" in st["not_implied"]
    c["revision_disclosed"] = "Q08-R1" in st["disclosed_revisions"]
    q7 = jr("certificates/Q07_P6_status.json")
    c["Q07_result_unchanged"] = q7["sigma0_closure_kappa3"].startswith("CLOSED at pole level")
    o1 = jr("certificates/Q03_O1_record.json")
    c["O1_record_unchanged"] = o1["O1_closure"] == "CLOSED (scoped)" and o1["O1_closure_confirmed"] is True and \
        st["O1_closure_unchanged"] == "CLOSED (scoped)" and st["bulk_symmetry_restored"] is False
    q7r = jr("certificates/Q07_release.json") if (CERT / "Q07_release.json").is_file() else {}
    c["Q07_layer_accepted"] = q7r.get("decision") == "ACCEPT"
    c["gate_reads_in_source_set"] = READS <= set(SS.source_files(Q)) | {f"certificates/{n}.json" for n in certs} | \
        {"certificates/Q03_O1_record.json", "certificates/Q07_release.json"}
    return c, st


def main():
    for p in STATUS:
        p.unlink(missing_ok=True)
    try:
        c, st = checks()
    except Exception as e:  # noqa: BLE001
        c, st = {"exception": False, "error": repr(e)}, {}
    ok = all(v is True for k, v in c.items() if k != "error")
    rec = {"landmark": "XVI-Q", "release": "Q08", "decision": "ACCEPT" if ok else "REJECT", "checks": c}
    (CERT / "Q08_release.json").write_text(json.dumps(rec, indent=1) + "\n")
    if not ok:
        print(json.dumps(rec, indent=1))
        return 1
    status = {"landmark": "XVI-Q", "release": "Q08", "status": "released (honest partial: P.6 evanescent CT structure)",
              "P6_outcome": st["P6_outcome"], "sigma0_closure_kappa3": st["sigma0_closure_kappa3"],
              "O1_closure": "CLOSED (scoped)", "O1_closure_confirmed": True, "hard_stops": HARD_STOPS,
              "independent_rebuild_q05": "not requested", "independent_rebuild_q06": "not requested",
              "independent_rebuild_q07": "in progress (reviewer)", "independent_rebuild_q08": "pending", "goals": st["goals"], "execution_pass": True, "claim_status": "partial", "gate": SELF}
    for p in STATUS:
        p.write_text(json.dumps(status, indent=1) + "\n")
    (Q / "ledgers").mkdir(exist_ok=True)
    (Q / "ledgers" / "Q08_milestones.json").write_text(json.dumps({"landmark": "XVI-Q", "release": "Q08", **strip(st)},
                                                                  indent=1) + "\n")
    print("ACCEPT")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
