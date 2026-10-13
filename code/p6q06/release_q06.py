#!/usr/bin/env python3
"""XVI-Q Q06 release gate (fails closed; on reject no RELEASE_STATUS.json is left).
(1) code/p6q06/q06_expected.json pinned; EVERY Q01..Q06 certificate equals its pinned copy; EVERY code/**/*.py,
    tools/*.py, tests/*.py file (incl. the Q05 gate) is pinned and re-hashed; an unpinned file is a reject;
(2) hardened static scan (p6q05/guard_scan2) run live over code/ (incl. p6q06): clean, reviewer evasions flagged;
(3) semantic: ghost rows closed on the 4D part with all mutations detected and a unique unconstrained fit, evanescent
    remainder explicit (every term has a T) and NOT claimed closed; G10 second route equal; V4 gauge identity exact;
    fast path equals Expr path; G01/G02 nonzero with Bose checks; h.h.c Missing; kappa^3 sigma0 closure NOT
    established; P6 Partial; O.1 record unchanged; Q05 layer gate accepted."""
from __future__ import annotations
import hashlib, json, sys  # noqa: E401
from pathlib import Path

Q = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Q / "tools"))
import source_set_q06 as SS  # noqa: E402
CERT = Q / "certificates"
EXPECTED_SHA256 = "1e8394983c4364b7192baf360ced8a93c1083229d07fc5764ba325f78f4da291"
STATUS = (Q / "RELEASE_STATUS.json", CERT / "RELEASE_STATUS.json")
HARD_STOPS = {"bulk_symmetry_restored": False, "insertion": "frozen", "XV_B_opened": False, "XIV_J_opened": False,
              "Route_B_discharged": False, "Missing_never_zero": True}
READS: set[str] = set()
Q06_CERTS = ["Q06_ghost_rows", "Q06_G10_crosscheck", "Q06_V4_validation", "Q06_fastpath", "Q06_P6_hhh", "Q06_P6_status"]


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
    # Q07-R1 (disclosed scope revision): Q07 files live under code/p6q06/q07/ and tools/*q06_q07*; they are pinned by
    # the Q07 gate (code/p6q06/q07/release_q07.py), not by this layer
    out = {f for f in out if not (f.startswith("code/p6q06/q07/") or (f.startswith("tools/") and "q07" in f))}
    return out - {"code/p6q06/release_q06.py"}


def checks():
    c = {}
    eb = rd("code/p6q06/q06_expected.json") if (Q / "code/p6q06/q06_expected.json").is_file() else b""
    c["expected_file_pinned"] = bool(eb) and sha_b(eb) == EXPECTED_SHA256
    exp = json.loads(eb) if eb else {}
    certs = exp.get("certificates", {})
    def _layer(n, o):
        # Q07-R1: the recursive scan count grows with later layers; the pinned count is a lower bound
        if n == "Q05_guard_scan" and isinstance(o, dict):
            o = dict(o)
            ok_n = isinstance(o.get("n_files_scanned"), int) and o["n_files_scanned"] >= certs[n]["n_files_scanned"]
            o["n_files_scanned"] = certs[n]["n_files_scanned"] if ok_n else -1
        return o
    bad = [n for n, ref in certs.items() if not (CERT / f"{n}.json").is_file()
           or _layer(n, strip(jr(f"certificates/{n}.json"))) != ref]
    c["every_certified_value_matches_pinned"] = len(certs) >= 30 and set(Q06_CERTS) <= set(certs) and not bad
    fs = exp.get("file_sha256", {})
    c["pinned_files_rehashed"] = bool(fs) and all((Q / f).is_file() and sha_b(rd(f)) == h for f, h in fs.items())
    c["no_unpinned_code_tools_tests_file"] = bool(fs) and code_files() <= set(fs)
    sys.path.insert(0, str(Q / "code" / "p6q05"))
    import guard_scan2
    READS.add("code/p6q05/guard_scan2.py")
    live = guard_scan2.report()
    c["scan_live_clean_incl_q06"] = live["scan_clean"] is True and live["all_reviewer_evasions_flagged"] is True and \
        "recursive" in live["scope"]
    r = jr("certificates/Q06_ghost_rows.json")
    ok = True
    for t in ("c.c.K", "c.c.cbar"):
        x = r[t]
        sol = x["unconstrained_fit_4D"]["solution_up_to_scale"]
        ok &= x["residual_4D_zero"] is True and x["residual_4D"] == "0" and len(x["mutations_detected"]) == 8 and \
            all(v is True for v in x["mutations_detected"].values()) and sol is not None and "0" not in sol and \
            x["evanescent_remainder_every_term_has_T"] is True and x["evanescent_remainder_n_terms"] > 0 and \
            x["evanescent_remainder_sha256"] == sha_b(x["evanescent_remainder"].encode()) and \
            x["channel_coefficients"]["K1"] == "-1" and x["channel_coefficients"]["K2"] == "1"
    c["ghost_rows_4D_closed_with_mutations"] = ok
    c["ghost_rows_evanescent_not_overclaimed"] = "OPEN" in r["status"] and "Grassmann" in r["diagnosis"]
    g10 = jr("certificates/Q06_G10_crosscheck.json")
    c["G10_second_route_equal"] = g10["route2_equals_Q05_G10"] is True
    v = jr("certificates/Q06_V4_validation.json")
    c["V4_gauge_identity_exact"] = v["gauge_n2_V3_vs_V2_zero"] is True and v["gauge_n3_V4_vs_V3_zero"] is True and \
        v["V4_permutation_symmetric"] is True
    f = jr("certificates/Q06_fastpath.json")
    c["fastpath_equals_expr_path"] = f == {"G06_fast_equals_expr": True, "G05_fast_equals_expr": True}
    h = jr("certificates/Q06_P6_hhh.json")["graphs"]
    c["G01_G02_evaluated"] = set(h) == {"G01", "G02"} and all(
        x["zero"] is False and x["sha256"] == sha_b(x["pole_over_U_d4_residue"].encode()) and x["bose_12_symmetric"] is True and x["bose_13_symmetric"] is True
        for x in h.values()) and h["G02"]["channels_bose_related"] is True
    st = jr("certificates/Q06_P6_status.json")
    c["honest_partial"] = st["P6_outcome"] == "Partial" and st["sigma0_closure_kappa3"].startswith("NOT") and \
        st["stages"]["graphs_missing"] == [] and st["stages"]["row_h_h_c"].startswith("Missing") and \
        st["kappa3_defect"] == "Missing" and st["kappa3_CT_solve"] == "Missing" and st["finite_parts"].startswith("Missing")
    c["errata_recorded"] = set(st["errata"]) == {"Q06-E1", "Q06-E2"} and "Q06-R1" in st.get("disclosed_revisions", {})
    o1 = jr("certificates/Q03_O1_record.json")
    c["O1_record_unchanged"] = o1["O1_closure"] == "CLOSED (scoped)" and o1["O1_closure_confirmed"] is True and \
        st["O1_closure_unchanged"] == "CLOSED (scoped)" and st["bulk_symmetry_restored"] is False
    q5 = jr("certificates/Q05_release.json") if (CERT / "Q05_release.json").is_file() else {}
    c["Q05_layer_accepted"] = q5.get("decision") == "ACCEPT"
    c["gate_reads_in_source_set"] = READS <= set(SS.source_files(Q)) | {f"certificates/{n}.json" for n in certs} | \
        {f"certificates/{n}.json" for n in Q06_CERTS} | {"certificates/Q03_O1_record.json", "certificates/Q05_release.json"}
    return c, st


def main():
    for p in STATUS:
        p.unlink(missing_ok=True)
    try:
        c, st = checks()
    except Exception as e:  # noqa: BLE001
        c, st = {"exception": False, "error": repr(e)}, {}
    ok = all(v is True for k, v in c.items() if k != "error")
    rec = {"landmark": "XVI-Q", "release": "Q06", "decision": "ACCEPT" if ok else "REJECT", "checks": c}
    (CERT / "Q06_release.json").write_text(json.dumps(rec, indent=1) + "\n")
    if not ok:
        print(json.dumps(rec, indent=1))
        return 1
    status = {"landmark": "XVI-Q", "release": "Q06", "status": "released (honest partial: P.6 continued)",
              "P6_outcome": st["P6_outcome"], "stages": st["stages"], "O1_closure": "CLOSED (scoped)",
              "O1_closure_confirmed": True, "hard_stops": HARD_STOPS,
              "independent_rebuild_passed_q04": True, "independent_rebuild_q05": "not requested (user: Q06 only)",
              "independent_rebuild_q06": "pending",
              "execution_pass": True, "claim_status": "partial", "gate": "code/p6q06/release_q06.py"}
    for p in STATUS:
        p.write_text(json.dumps(status, indent=1) + "\n")
    (Q / "ledgers").mkdir(exist_ok=True)
    (Q / "ledgers" / "Q06_milestones.json").write_text(json.dumps({"landmark": "XVI-Q", "release": "Q06",
        "P6_outcome": st["P6_outcome"], "stages": st["stages"], "hard_stops": HARD_STOPS}, indent=1) + "\n")
    print("ACCEPT")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
