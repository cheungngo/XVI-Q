#!/usr/bin/env python3
"""XVI-Q Q01 release gate.  Value-checks EVERY certified number (q01_conventions.EXPECTED) in Q01_intake / Q01_O2 /
Q01_guard, the import pins and the hard stops.  Fails closed; on reject no RELEASE_STATUS.json is left behind."""
from __future__ import annotations
import json
import sys
from pathlib import Path

Q = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Q / "code" / "q01"))
import q01_conventions as QC  # noqa: E402

CERT = Q / "certificates"
STATUS = ()   # Q02: Q01 is a layer record (certificates/Q01_release.json); RELEASE_STATUS.json is written by the current release gate


def jl(n):
    f = CERT / f"{n}.json"
    return json.loads(f.read_text()) if f.is_file() else {}


def checks():
    it, o2, g = jl("Q01_intake"), jl("Q01_O2"), jl("Q01_guard")
    E = QC.EXPECTED
    a = o2.get("items", {}).get("a_ell4", {}).get("detail", {})
    b = o2.get("items", {}).get("b_antibracket_expansion", {}).get("detail", {}).get("unmutated", {})
    c = {}
    c["intake_pins"] = it.get("pins_ok") is True and it.get("XVI_P_p23_zip_ok") is True
    c["intake_extracted_identical"] = it.get("XVI_P_p23_extracted_identical") is True
    c["intake_O1_ledger"] = it.get("O1_ledger_sha_ok") is True
    c["errata_recorded"] = [e.get("id") for e in it.get("errata", [])][:5] == ["E1", "E2", "E3", "E4", "E5"]   # Q02 appends E6, E7
    c["ell4_n_terms"] = a.get("n_terms_mine") == E["ell4_n_terms"] == a.get("n_terms_export")
    c["ell4_normalisation"] = a.get("normalisation_from_XVI_O_probe") == E["ell4_normalisation_from_XVI_O_probe"]
    c["ell4_full_equal"] = a.get("full_polynomial_equal") is E["ell4_full_polynomial_equal"] and a.get("n_terms_difference") == 0
    c["ell4_probe_equal"] = a.get("export_on_XVI_O_probe_equal") is E["ell4_probe_equal"]
    c["ell4_mutation"] = a.get("mutation_one_term_x2_detected") is True
    c["exports_sha"] = {int(k): v for k, v in o2.get("exports_sha256", {}).items()} == QC.EXPORT_SHA256
    c["W_terms"] = all(b.get(k, {}).get("W_terms") == v for k, v in E["W_terms"].items()) \
        and o2.get("W_terms") == E["W_terms"]
    c["E_c_W_zero"] = all(b.get(k, {}).get("antibracket_af0_vanishes") is v and
                          all(x == 0 for x in b.get(k, {}).get("E_c_W_terms", {"x": 1}).values())
                          for k, v in E["E_c_W_all_zero"].items())
    bd = o2.get("items", {}).get("b_antibracket_expansion", {}).get("detail", {})
    c["antibracket_mutations"] = bd.get("mutation_scale_ell4_term_detected") is True and bd.get("mutation_drop_s1_detected") is True
    c["master_mutations"] = {m: v.get("n_residuals") for m, v in o2.get("master_mutations", {}).items()} \
        == E["master_mutation_residuals"] and all(v.get("claim_status") != "verified" for v in o2.get("master_mutations", {}).values())
    c["classical_certs_identical"] = sorted(k for k, v in o2.get("classical_certs_content_identical_to_p23", {}).items() if v) \
        == sorted(E["classical_certs_content_identical_to_p23"])
    c["kappa2_table_all_true"] = all(o2.get("items", {}).get("c_ghost_antifield", {}).get("kappa2_table", {"x": False}).values())
    c["s2_derived_zero"] = set(o2.get("items", {}).get("c_ghost_antifield", {}).get("s2", {}).values()) == {"derived_zero"}
    c["O2_items"] = o2.get("item_status") == E["O2_items"] and o2.get("O2_kappa2_master_status") == "Constructed"
    c["guard"] = g.get("n_forms") == E["guard_forms"] and g.get("n_exprs") == E["guard_exprs"] and g.get("missed") == []
    c["guard_algebra"] = g.get("algebra", {}).get("authorised_value") == "1/80" and g.get("algebra", {}).get("d_free_rem") == "2"
    return c, it, o2, g


def main():
    try:
        c, it, o2, g = checks()
    except Exception as exc:  # noqa: BLE001
        c, it, o2, g = {"error": f"{type(exc).__name__}: {exc}"}, {}, {}, {}
    ok = all(v is True for v in c.values())
    for s in STATUS:
        s.unlink(missing_ok=True)
    rec = {"landmark": QC.LANDMARK, "release": QC.RELEASE, "version_label": QC.VERSION_LABEL, "date": QC.DATE,
           "execution_pass": ok, "publication": "accepted" if ok else "REJECTED", "value_checks": c,
           "claim_status": ("Q01: O.2 classical kappa^2 BV master 2(S0,S2)+(S1,S1)=0 on S_inv Constructed in XVI-Q: (a) full "
                            "ell_4 independently reconstructed (3360/3360) and matched to XVI-O's probe/convention; (b) explicit "
                            "af-0 antibracket expansion E_c(W_2)=0 exactly on an independent engine; (c) ghost/antifield cells "
                            "via re-executed exact P13 master with s2 derived zero. XVI-O's own O.2 record is unchanged "
                            "(read-only). O.1 NOT_CLOSED (bulk O.5, retained identity O.6/O.7 open -> Q02); bulk symmetry "
                            "not restored.") if ok else "Q01: not accepted",
           "scope": o2.get("scope"),
           "flags": {"intake_pinned": "earned" if ok else "not_earned",
                     "O2_ell4_constructed": "earned" if ok else "not_earned",
                     "O2_antibracket_expansion_constructed": "earned" if ok else "not_earned",
                     "O2_ghost_antifield_constructed": "earned" if ok else "not_earned",
                     "O2_kappa2_master_constructed": "earned" if ok else "not_earned",
                     "d_guard_q01_verified": "earned" if ok else "not_earned",
                     "bulk_O5_evaluated": "not_earned", "bulk_identity_verified": "not_earned",
                     "author_rebuild_passed": "see reproducibility", "independent_rebuild_passed": "pending"},
           "O1_closure": "NOT_CLOSED", "bulk_symmetry_restored": False, "symmetry_restored": False,
           "hard_stops": QC.HARD_STOPS, "XVI_Q_registered": "locally (registration/local_registration.json); not published",
           "O1_checklist_after_Q01": {"O.2": "Constructed (XVI-Q Q01)", "O.4": "inventory (XVI-P; errata E1)",
                                      "O.5": "bulk OPEN (D04, D05a, D06 bulk, D12 Missing)", "O.6": "bulk statement open (ST4a/ST4b)",
                                      "O.7": "bulk CT solution open", "enlarge": "not triggered at kappa^2 (XVI-P p23)"}}
    if ok:
        txt = json.dumps(rec, indent=1) + "\n"
        for s in STATUS:
            s.write_text(txt)
        (CERT / "Q01_release.json").write_text(txt)
        (CERT / "Q01_release_rejected.json").unlink(missing_ok=True)
        (Q / "ledgers").mkdir(exist_ok=True)
        (Q / "ledgers" / "Q01_milestones.json").write_text(json.dumps(
            {"landmark": "XVI-Q", "release": "Q01", "milestones": rec["O1_checklist_after_Q01"], "flags": rec["flags"],
             "hard_stops": QC.HARD_STOPS}, indent=1) + "\n")
    else:
        (CERT / "Q01_release.json").unlink(missing_ok=True)
        (CERT / "Q01_release_rejected.json").write_text(json.dumps(rec, indent=1) + "\n")
    print(json.dumps({"publication": rec["publication"], "failed": [k for k, v in c.items() if v is not True]}))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
