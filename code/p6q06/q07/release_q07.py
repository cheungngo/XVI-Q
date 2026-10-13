#!/usr/bin/env python3
"""XVI-Q Q07 release gate (fails closed; on reject no RELEASE_STATUS.json is left).
(1) code/p6q06/q07/q07_expected.json pinned; EVERY Q01..Q07 certificate equals its pinned copy (Q05 guard-scan
    file count: lower bound); EVERY code/**/*.py, tools/*.py, tests/*.py (incl. the Q05 and Q06 gates) pinned and
    re-hashed; an unpinned file is a reject;
(2) live hardened static scan clean;
(3) semantic: h.h.c tree calibrated, 4D residual 0, every mutation detected, unique fit, matter sector exact,
    remainder explicit and evanescent; evanescent-CT naive test recorded as NOT consistent; joint local CT solvable,
    verified, Bose, evanescent, CT shas; status consistent; O.1 unchanged; Q06 layer accepted; Q07-R1 disclosed."""
from __future__ import annotations
import hashlib, json, sys  # noqa: E401
from pathlib import Path

Q = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Q / "tools"))
import source_set_q06_q07 as SS  # noqa: E402
CERT = Q / "certificates"
EXPECTED_SHA256 = "8d2c4d5c63816e89e1b883a2cd63df6d1c3457ebe914790c2594a68ff980e746"
STATUS = (Q / "RELEASE_STATUS.json", CERT / "RELEASE_STATUS.json")
HARD_STOPS = {"bulk_symmetry_restored": False, "insertion": "frozen", "XV_B_opened": False, "XIV_J_opened": False,
              "Route_B_discharged": False, "Missing_never_zero": True}
READS: set[str] = set()
Q07_CERTS = ["Q07_hhc_row", "Q07_evanescent_CT", "Q07_P6_status"]
SELF = "code/p6q06/q07/release_q07.py"


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
    # Q08-R1 (disclosed scope revision): Q08 files (code/p6q06/q07/q08/, tools/*q08*) are pinned by the Q08 gate
    out = {f for f in out if not (f.startswith("code/p6q06/q07/q08/") or (f.startswith("tools/") and "q08" in f))}
    return out - {SELF}


def checks():
    c = {}
    ep = "code/p6q06/q07/q07_expected.json"
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
    c["every_certified_value_matches_pinned"] = len(certs) >= 35 and set(Q07_CERTS) <= set(certs) and not bad
    fs = exp.get("file_sha256", {})
    c["pinned_files_rehashed"] = bool(fs) and all((Q / f).is_file() and sha_b(rd(f)) == h for f, h in fs.items())
    c["no_unpinned_code_tools_tests_file"] = bool(fs) and code_files() <= set(fs) and \
        {"code/p6q05/release_q05.py", "code/p6q06/release_q06.py"} <= set(fs)
    sys.path.insert(0, str(Q / "code" / "p6q05"))
    import guard_scan2
    READS.add("code/p6q05/guard_scan2.py")
    live = guard_scan2.report()
    c["scan_live_clean_incl_q07"] = live["scan_clean"] is True and live["all_reviewer_evasions_flagged"] is True
    h = jr("certificates/Q07_hhc_row.json")
    c["hhc_tree_calibrated"] = h["tree_calibration"] == {"sigma_F": "I", "tree_row_closes": True,
                                                         "tree_row_without_E_nonzero": True}
    sol = h["unconstrained_fit_4D"]["solution"]
    c["hhc_4D_closed_with_mutations"] = h["residual_4D_zero"] is True and h["residual_4D"] == "0" and \
        len(h["mutations_detected"]) == 6 and all(v is True for v in h["mutations_detected"].values()) and \
        sol is not None and len(set(sol)) == 1 and sol[0] != "0" and h["matter_sector_closes_exactly"] is True
    c["hhc_remainder_explicit_evanescent"] = h["evanescent_remainder_every_term_has_T"] is True and \
        h["evanescent_remainder_n_terms"] > 0 and h["evanescent_remainder_sha256"] == sha_b(h["evanescent_remainder"].encode())
    e = jr("certificates/Q07_evanescent_CT.json")
    nv = e["naive_p22_type"]
    c["naive_CT_recorded_not_consistent"] = nv["ccK"]["zero"] is False and nv["ccb"]["zero"] is False and \
        nv["verdict"].startswith("NOT")
    j = e["joint_local_CT"]
    c["joint_CT_solvable_verified"] = j["solvable"] is True and j["verified_rows_zero"] is True and \
        j["hhh_CT_bose_symmetric"] is True and j["all_CT_evanescent"] is True and j["n_free"] > 0 and \
        all(j[f"particular_CT_{n}_sha256"] == sha_b(j[f"particular_CT_{n}"].encode())
            for n in ("Gamma_Kch", "Gamma_cbar_ch", "Gamma_cstar_cc", "Gamma_hhh")) and \
        e["verdict"].startswith("CANCELLABLE") and "not unique" in e["verdict"]
    st = jr("certificates/Q07_P6_status.json")
    c["status_scoped_honest"] = st["P6_outcome"] == "Partial" and st["flags"]["claim_status"] == "partial" and \
        st["sigma0_closure_kappa3"].startswith("CLOSED at pole level on the kappa^3 one-loop THREE-POINT rows") and \
        "kappa^3 four-point / higher rows" in st["not_implied"] and "finite parts" in st["not_implied"] and \
        st["kappa3_CT_domain_reduction"].startswith("Partial")
    c["revision_disclosed"] = "Q07-R1" in st["disclosed_revisions"] and "Q06-E2" in st["errata"]
    o1 = jr("certificates/Q03_O1_record.json")
    c["O1_record_unchanged"] = o1["O1_closure"] == "CLOSED (scoped)" and o1["O1_closure_confirmed"] is True and \
        st["O1_closure_unchanged"] == "CLOSED (scoped)" and st["bulk_symmetry_restored"] is False
    q6 = jr("certificates/Q06_release.json") if (CERT / "Q06_release.json").is_file() else {}
    c["Q06_layer_accepted"] = q6.get("decision") == "ACCEPT"
    c["gate_reads_in_source_set"] = READS <= set(SS.source_files(Q)) | {f"certificates/{n}.json" for n in certs} | \
        {"certificates/Q03_O1_record.json", "certificates/Q06_release.json"}
    return c, st


def main():
    for p in STATUS:
        p.unlink(missing_ok=True)
    try:
        c, st = checks()
    except Exception as e:  # noqa: BLE001
        c, st = {"exception": False, "error": repr(e)}, {}
    ok = all(v is True for k, v in c.items() if k != "error")
    rec = {"landmark": "XVI-Q", "release": "Q07", "decision": "ACCEPT" if ok else "REJECT", "checks": c}
    (CERT / "Q07_release.json").write_text(json.dumps(rec, indent=1) + "\n")
    if not ok:
        print(json.dumps(rec, indent=1))
        return 1
    status = {"landmark": "XVI-Q", "release": "Q07", "status": "released (honest partial: P.6 three-point rows)",
              "P6_outcome": st["P6_outcome"], "stages": st["stages"], "sigma0_closure_kappa3": st["sigma0_closure_kappa3"],
              "O1_closure": "CLOSED (scoped)", "O1_closure_confirmed": True, "hard_stops": HARD_STOPS,
              "independent_rebuild_q05": "not requested", "independent_rebuild_q06": "not requested",
              "independent_rebuild_q07": "pending", "execution_pass": True, "claim_status": "partial", "gate": SELF}
    for p in STATUS:
        p.write_text(json.dumps(status, indent=1) + "\n")
    (Q / "ledgers").mkdir(exist_ok=True)
    (Q / "ledgers" / "Q07_milestones.json").write_text(json.dumps({"landmark": "XVI-Q", "release": "Q07", **strip(st)},
                                                                  indent=1) + "\n")
    print("ACCEPT")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
