#!/usr/bin/env python3
"""XVI-Q Q03 release gate (fails closed; on reject no RELEASE_STATUS.json is left).
(1) code/q03/q03_expected.json is pinned (EXPECTED_SHA256) and EVERY Q01/Q02/Q03 certificate equals its pinned copy;
(2) RE-HASH at gate time: every pinned provenance file in intake/imports.json (no unpinned provenance file allowed),
    the bound Q02 independent report, and every pinned code/conventions file;
(3) every file the gate reads is either in the source set (tools/source_set.py) or a build output it re-checks;
(4) semantic checks of the O.1 record (enlarged domain, prescription, reason 8 not answered, binding, scope, not-implied,
    caveats) and of the ring/field guard forms; the Q01 and Q02 layer gates accepted."""
from __future__ import annotations
import hashlib, json, sys  # noqa: E401
from pathlib import Path

Q = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Q / "tools"))
import source_set as SS  # noqa: E402
CERT = Q / "certificates"
EXPECTED_FILE = Q / "code" / "q03" / "q03_expected.json"
EXPECTED_SHA256 = "44cdd33648bd8ba42267fbd70f28a15d6c32dd988b64be8170204760e23aa9b0"
STATUS = (Q / "RELEASE_STATUS.json", CERT / "RELEASE_STATUS.json")
if "--layer" in sys.argv:   # Q04: earlier-layer gate; RELEASE_STATUS is owned by the Q04 gate
    STATUS = ()
GENERATED = ("certificates/", "ledgers/", "build/")       # build outputs: checked by value / re-hash, not source
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


def checks():
    c = {}
    eb = rd("code/q03/q03_expected.json") if EXPECTED_FILE.is_file() else b""
    c["expected_file_pinned"] = bool(eb) and sha_b(eb) == EXPECTED_SHA256
    exp = json.loads(eb) if eb else {}
    bad = []
    for name, ref in exp.get("certificates", {}).items():
        f = f"certificates/{name}.json"
        if not (Q / f).is_file() or strip(jr(f)) != ref:
            bad.append(name)
    c["every_certified_value_matches_pinned"] = len(exp.get("certificates", {})) >= 12 and not bad
    c["_value_mismatches"] = bad
    badf = [f for f, s in exp.get("file_sha256", {}).items() if not (Q / f).is_file() or sha_b(rd(f)) != s]
    c["code_and_conventions_rehashed"] = len(exp.get("file_sha256", {})) >= 8 and not badf
    c["_file_mismatches"] = badf
    imp = jr("intake/imports.json")
    badp = [f for f, m in imp["files"].items() if not (Q / f).is_file() or sha_b(rd(f)) != m["sha256"]]
    unpinned = sorted(p.relative_to(Q).as_posix() for p in (Q / "provenance").rglob("*") if p.is_file()
                      and p.relative_to(Q).as_posix() not in imp["files"])
    c["provenance_rehashed_all_pinned"] = len(imp["files"]) == 33 and not badp and not unpinned
    c["_provenance_mismatches"] = badp + [f"unpinned:{u}" for u in unpinned]
    rec = jr("certificates/Q03_O1_record.json")
    ind = "reproducibility/independent_rebuild_report_q02.json"
    rb = rd(ind)
    cf = rec["confirmation"]
    c["bound_independent_report_rehashed"] = sha_b(rb) == cf["report_sha256"] == "e4b54702bc768899f22a00a289de0c0aeb2bf9317b87a28c4a462f05ae4bf495" \
        and json.loads(rb)["O1_verdict"]["decision"] == "CONFIRM CLOSED (scoped)" and cf["binding_ok"] is True \
        and json.loads(rb)["source"]["tree_sha256"] == cf["reviewed_tree_sha256"]
    dom = rec["counterterm_domain"]
    c["item5_enlarged_domain"] = dom["n_p14r_4D"] == 11 and dom["n_p22_evanescent"] == 10 and dom["n_enlarged_domain"] == 21 \
        and dom["p22_evanescent_coefficients"]["Kc_E2"] == "-1/16" and dom["p22_evanescent_coefficients"]["ccb_SE"] == "1/8" \
        and any("ENLARGED domain" in i["evidence"] and "p22 evanescent" in i["evidence"] for i in rec["items"]
                if i["item"] == "Demonstrate identity remainder cancellation / CT solution")
    r = rec["closure_reasons"]
    r8 = list(r.values())[7]
    c["reason8_not_answered_prescription_stated"] = len(r) == 8 and not r8["status"].startswith("ANSWERED") \
        and "XVI-P" in rec["prescription"]["used"] and "o70-sync" in rec["prescription"]["not_used"] \
        and all(v["status"] == "ANSWERED" for v in list(r.values())[:7])
    c["scope_verbatim"] = rec["scope_verbatim_from_independent_review"] == json.loads(rb)["O1_verdict"]["scope"]
    ni = " | ".join(rec["not_implied"])
    c["not_implied_complete"] = all(s in ni for s in ("ST4b", "bulk symmetry restoration", "XVI-O", "NOT_CLOSED", "P.6", "insertion"))
    cv = " | ".join(rec["caveats"])
    c["caveats_recorded"] = all(s in cv for s in ("same-box", "D01/D03", "sigma_0 closure", "c.h defect"))
    c["O1_closed_scoped_consistent"] = rec["O1_closure"] == "CLOSED (scoped)" and rec["O1_closure_confirmed"] is True \
        and rec["all_earned"] is True and len(rec["items"]) == 6 and all(i["earned"] is True for i in rec["items"])
    g = jr("certificates/Q03_guard_ring.json")
    c["guard_ring_field_forms_all_raise"] = g["n_forms"] == 8 and g["missed"] == []
    g1 = jr("certificates/Q01_guard.json")
    c["guard_18_forms_still_raise"] = g1["n_forms"] == 18 and g1["missed"] == []
    c["Q01_layer_accepted"] = jr("certificates/Q01_release.json").get("publication") == "accepted"
    q2 = jr("certificates/Q02_release.json")
    c["Q02_layer_accepted"] = q2.get("publication") == "accepted" and q2.get("bulk_symmetry_restored") is False
    for f in ("intake/errata_p23.json", "registration/local_registration.json", "code/q02/q02_expected.json"):
        rd(f)
    src = set(SS.source_files(Q))
    outside = sorted(f for f in READS if f not in src and not f.startswith(GENERATED))
    c["every_gate_read_in_source_set"] = not outside and len(READS) > 40
    c["_reads_outside_source_set"] = outside
    return c, rec


def main():
    try:
        c, rec = checks()
    except Exception as exc:  # noqa: BLE001
        c, rec = {"error": f"{type(exc).__name__}: {exc}"}, {}
    ok = all(v is True for k, v in c.items() if not k.startswith("_"))
    for s in STATUS:
        s.unlink(missing_ok=True)
    closure = rec.get("O1_closure", "NOT_CLOSED") if ok else "NOT_CLOSED"
    out = {"landmark": "XVI-Q", "release": "Q03", "version_label": "XVI-Q Q03", "date": "2026-10-10",
           "kind": "correction release (no new science)", "execution_pass": ok,
           "publication": "accepted" if ok else "REJECTED", "value_checks": c,
           "O1_closure": closure, "O1_closure_confirmed": bool(ok and rec.get("O1_closure_confirmed")),
           "O1_confirmation": rec.get("confirmation"), "O1_closure_scope": rec.get("scope_verbatim_from_independent_review"),
           "O1_prescription": rec.get("prescription"), "O1_counterterm_domain": rec.get("counterterm_domain", {}).get("statement"),
           "O1_not_implied": rec.get("not_implied"), "O1_caveats": rec.get("caveats"),
           "O1_checklist": [{k: i[k] for k in ("item", "owner", "earned", "evidence")} for i in rec.get("items", [])],
           "XVI_O_ledger_O1": "NOT_CLOSED (read-only, unchanged)", "bulk_symmetry_restored": False, "symmetry_restored": False,
           "flags": {"O1_closed_scoped": "earned + independently confirmed (Q02 review)" if ok else "not_earned",
                     "O2_kappa2_master_constructed": "earned (Q01)", "bulk_identity_verified_two_point_kappa2": "earned (Q02)",
                     "d_guard_q03_ring_field": "earned" if ok else "not_earned",
                     "independent_rebuild_passed_q01": True, "independent_rebuild_passed_q02": True,
                     "independent_rebuild_passed_q03": "pending", "author_rebuild_passed": "see reproducibility"},
           "hard_stops": HARD_STOPS, "XVI_Q_registered": "locally; not published"}
    cr = CERT / "Q03_release.json"
    if ok:
        txt = json.dumps(out, indent=1) + "\n"
        for s in STATUS:
            s.write_text(txt)
        cr.write_text(txt)
        (CERT / "Q03_release_rejected.json").unlink(missing_ok=True)
        (Q / "ledgers" / "Q03_milestones.json").write_text(json.dumps(
            {"landmark": "XVI-Q", "release": "Q03", "O1_closure": closure, "O1_closure_confirmed": out["O1_closure_confirmed"],
             "flags": out["flags"], "hard_stops": HARD_STOPS}, indent=1) + "\n")
    else:
        cr.unlink(missing_ok=True)
        (CERT / "Q03_release_rejected.json").write_text(json.dumps(out, indent=1, default=str) + "\n")
    print(json.dumps({"publication": out["publication"], "O1_closure": closure,
                      "failed": [k for k, v in c.items() if not k.startswith("_") and v is not True],
                      "mm": {k: v for k, v in c.items() if k.startswith("_") and v}}, default=str)[:1500])
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
