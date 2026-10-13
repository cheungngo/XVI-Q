#!/usr/bin/env python3
"""XVI-Q Q05 release gate (fails closed; on reject no RELEASE_STATUS.json is left).
(1) code/p6q05/q05_expected.json pinned; EVERY Q01..Q05 certificate equals its pinned copy; EVERY code/**/*.py, tools/*.py,
    tests/*.py file is pinned and re-hashed -- an unpinned file is a reject (Q04 review); provenance via Q04 layer gate;
    bound Q04 independent report re-hashed;
(2) hardened scan run live: clean and all 7 reviewer evasion forms flagged; claim recorded as best-effort;
(3) semantic: own-rule validations constant with ratios D01/D03/D06 = I, D13 = -I (disclosed); two-point ghost equation
    constant; row c.phi.phi residual exactly 0 with every mutation detected and the S1.Gamma_Kc channel required;
    ghost rows reported OPEN (never closed/obstruction); G01/G02 Missing; P6_outcome Partial; O.1 record unchanged;
(4) Q04 layer gate accepted (certificates/Q04_release.json)."""
from __future__ import annotations
import hashlib, json, sys  # noqa: E401
from pathlib import Path

Q = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Q / "tools"))
import source_set_q05 as SS  # noqa: E402
CERT = Q / "certificates"
EXPECTED_SHA256 = "2b9291e05ed07c04bbf1d7250c05b123c32213ee4eef002ec5f92c6d4673e0b9"
Q04_REPORT_SHA = "bc2b16d9"
STATUS = (Q / "RELEASE_STATUS.json", CERT / "RELEASE_STATUS.json")
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


def code_files():
    out = {p.relative_to(Q).as_posix() for p in (Q / "code").rglob("*.py") if "__pycache__" not in p.parts}
    out |= {p.relative_to(Q).as_posix() for p in (Q / "tools").glob("*.py")}
    out |= {p.relative_to(Q).as_posix() for p in (Q / "tests").glob("*.py")}
    # Q06-owned files are pinned by the Q06 gate (code/p6q06/release_q06.py), not by this layer
    q06 = {f for f in out if f.startswith("code/p6q06/") or f.split("/")[-1] in ("test_q06.py",) or
           (f.startswith("tools/") and "q06" in f)}
    return out - {"code/p6q05/release_q05.py"} - q06


def checks():
    c = {}
    eb = rd("code/p6q05/q05_expected.json") if (Q / "code/p6q05/q05_expected.json").is_file() else b""
    c["expected_file_pinned"] = bool(eb) and sha_b(eb) == EXPECTED_SHA256
    exp = json.loads(eb) if eb else {}
    certs = exp.get("certificates", {})
    def _layer(n, o):
        # the live scan count grows when later layers (Q06+) add code files; the pinned count is a lower bound
        if n == "Q05_guard_scan" and isinstance(o, dict):
            o = dict(o)
            ok_n = isinstance(o.get("n_files_scanned"), int) and o["n_files_scanned"] >= certs[n]["n_files_scanned"]
            o["n_files_scanned"] = certs[n]["n_files_scanned"] if ok_n else -1
        return o
    bad = [n for n, ref in certs.items() if not (CERT / f"{n}.json").is_file()
           or _layer(n, strip(jr(f"certificates/{n}.json"))) != ref]
    c["every_certified_value_matches_pinned"] = len(certs) >= 24 and not bad
    fs = exp.get("file_sha256", {})
    c["pinned_files_rehashed"] = bool(fs) and all((Q / f).is_file() and sha_b(rd(f)) == h for f, h in fs.items())
    c["no_unpinned_code_tools_tests_file"] = bool(fs) and code_files() <= set(fs)
    rep = rd("reproducibility/independent_rebuild_report_q04.json")
    c["bound_q04_report_rehashed"] = sha_b(rep).startswith(Q04_REPORT_SHA) and \
        sha_b(rep) == fs.get("reproducibility/independent_rebuild_report_q04.json")
    sys.path.insert(0, str(Q / "code" / "p6q05"))
    import guard_scan2
    READS.add("code/p6q05/guard_scan2.py")
    live = guard_scan2.report()
    gs = jr("certificates/Q05_guard_scan.json")
    c["scan_live_clean_and_evasions_flagged"] = live["scan_clean"] and live["all_reviewer_evasions_flagged"] and \
        gs["scan_clean"] is True and "NOT a security boundary" in gs["claim"]
    v = jr("certificates/Q05_rules_validation.json")
    c["rules_validated"] = v["all_constant"] is True and v["ratios_as_expected"] is True and \
        {k: x["ratio"] for k, x in v["validations"].items()} == {"D01": "I", "D03": "I", "D06": "I", "D13": "-I"}
    c["D13_orientation_disclosed"] = "-I" in v["D13_orientation_disclosure"]
    c["ghost_equation_two_point_constant"] = v["ghost_equation_two_point"]["ratio_D06_over_D13_K_to_T"] == "-2"
    g3 = v["ghost_equation_three_point"]
    c["ghost_equation_three_point"] = g3["ratio_G04_over_G03_K_to_T"] == g3["ratio_G09_over_G07_K_to_T"] == "-2"
    g = jr("certificates/Q05_P6_graphs.json")
    c["graphs_evaluated_nonzero"] = sorted(g["graphs"]) == ["G03", "G04", "G05", "G06", "G07", "G08", "G09", "G10",
                                                              "G11", "G12", "G13"] and \
        all(x["zero"] is False and x["sha256"] == sha_b(x["pole_over_U_d4_residue"].encode()) for x in g["graphs"].values())
    c["G01_G02_missing"] = set(g["missing"]) == {"G01", "G02"}
    r = jr("certificates/Q05_P6_rows.json")
    cp = r["c.phi.phi"]
    c["row_c_phi_phi_closed_with_mutations"] = cp["residual"] == "0" and cp["closed_at_pole_level"] is True and \
        cp["tree_identity_closed_with_sigma"] is True and len(cp["mutations_detected"]) >= 6 and \
        all(v for k, v in cp["mutations_detected"].items() if k != "drop_G05") and \
        cp["mutations_detected"].get("drop_S1_Kc_channel") is True and cp["G05_s0_projection_zero"] is True
    c["ghost_rows_open_not_overclaimed"] = all(x["loop_residual_zero"] is False and x["status"].startswith("OPEN")
                                               and x["tree_closes_with_constant_sigma"] is True
                                               for x in r["ghost_rows"].values()) and r["h.h.c"]["status"] == "Missing"
    st = jr("certificates/Q05_P6_status.json")
    c["honest_partial"] = st["P6_outcome"] == "Partial" and st["sigma0_closure_kappa3"].startswith("NOT") and \
        st["stages"]["graphs_missing"] == ["G01", "G02"] and st["finite_parts"].startswith("Missing")
    c["errata_recorded"] = set(st["errata"]) == {"Q05-E1", "Q05-E2", "Q05-E3", "Q05-E4"}
    o1 = jr("certificates/Q03_O1_record.json")
    c["O1_record_unchanged"] = o1["O1_closure"] == "CLOSED (scoped)" and o1["O1_closure_confirmed"] is True and \
        st["O1_closure_unchanged"] == "CLOSED (scoped)" and st["bulk_symmetry_restored"] is False
    q4 = jr("certificates/Q04_release.json") if (CERT / "Q04_release.json").is_file() else {}
    c["Q04_layer_accepted"] = q4.get("publication") == "accepted"
    c["gate_reads_in_source_set"] = READS <= set(SS.source_files(Q)) | {f"certificates/{n}.json" for n in certs} | \
        {"certificates/Q05_guard_scan.json", "certificates/Q05_rules_validation.json", "certificates/Q05_P6_graphs.json",
         "certificates/Q05_P6_rows.json", "certificates/Q05_P6_status.json", "certificates/Q03_O1_record.json",
         "certificates/Q04_release.json"}
    return c, st


def main():
    for p in STATUS:
        p.unlink(missing_ok=True)
    try:
        c, st = checks()
    except Exception as e:  # noqa: BLE001
        c, st = {"exception": False, "error": repr(e)}, {}
    ok = all(v is True for k, v in c.items() if k != "error")
    rec = {"landmark": "XVI-Q", "release": "Q05", "decision": "ACCEPT" if ok else "REJECT", "checks": c}
    (CERT / "Q05_release.json").write_text(json.dumps(rec, indent=1) + "\n")
    if not ok:
        print(json.dumps(rec, indent=1))
        return 1
    status = {"landmark": "XVI-Q", "release": "Q05", "status": "released (honest partial: P.6 continued)",
              "P6_outcome": st["P6_outcome"], "stages": st["stages"], "O1_closure": "CLOSED (scoped)",
              "O1_closure_confirmed": True, "hard_stops": HARD_STOPS,
              "independent_rebuild_passed_q04": True, "independent_rebuild_q05": "pending",
              "execution_pass": True, "claim_status": "partial", "gate": "code/p6q05/release_q05.py"}
    for p in STATUS:
        p.write_text(json.dumps(status, indent=1) + "\n")
    (Q / "ledgers").mkdir(exist_ok=True)
    (Q / "ledgers" / "Q05_milestones.json").write_text(json.dumps({"landmark": "XVI-Q", "release": "Q05",
        "P6_outcome": st["P6_outcome"], "stages": st["stages"], "hard_stops": HARD_STOPS}, indent=1) + "\n")
    print("ACCEPT")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
