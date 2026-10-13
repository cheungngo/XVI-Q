#!/usr/bin/env python3
"""XVI-Q Q03 correction layer (no new science).  Writes
  certificates/Q03_O1_record.json -- the O.1 closure record: CLOSED (scoped), confirmed, bound to the Q02 independent
      review (report sha256 + reviewed source zip/tree); item-5 evidence with the ENLARGED counterterm domain (p14r 11 4D
      directions + p22 evanescent operators); prescription statement (XVI-P DRED, not XVI-O's o70-sync D11 kernel
      prescription; closure reason 8 NOT marked answered); scope verbatim; what it does not imply; caveats.
  certificates/Q03_guard_ring.json -- the PolyRing/FracField-with-d-generator bypass forms (Q02 review), each must raise.
Q02 certificates are a frozen layer and are not rewritten."""
from __future__ import annotations
import hashlib, json, sys  # noqa: E401
from pathlib import Path

Q = Path(__file__).resolve().parents[2]
CERT = Q / "certificates"
IND = Q / "reproducibility" / "independent_rebuild_report_q02.json"
BOUND = {"report_sha256": "e4b54702bc768899f22a00a289de0c0aeb2bf9317b87a28c4a462f05ae4bf495",
         "reviewed_zip_sha256": "378651a3a279cde0a8fd67d4e39ac5833d5a94ff00f10f62a001be8ab2b244dc",
         "reviewed_tree_sha256": "4b551b9801cd14b8f2b81775942d535ab282a75c1557148a94a91a9dad88a810",
         "decision": "CONFIRM CLOSED (scoped)"}
NOT_IMPLIED = [
    "the functional / all-n-point ST4b identity (only the retained two-point projection is verified)",
    "bulk symmetry restoration (bulk_symmetry_restored stays false)",
    "any change to XVI-O's own O.1 ledger: XVI-O is read-only and its O.1 stays NOT_CLOSED there",
    "the three-point sector / P.6 (kappa^3) and higher orders or loops",
    "the insertion (frozen; XVI-N massless obstruction untouched)",
    "finite renormalised amplitudes with non-QdS (epsilon-scalar) external components (p22 Q4S-trace O(eps)x pole remnant stays disclosed)"]
CAVEATS = [
    "same-box reviews: Q01/Q02 independent rebuilds ran on the author's machine (shared Python/sympy pin and TeX Live); not an independent machine",
    "D01/D03 pole values were never recomputed from independent Feynman rules in any review; cross-checked only via the exact cancellation of the c.h defect in d",
    "full sigma_0 closure and the c.h defect were not re-derived in the Q02 review; they rely on the earlier p19-p22 independent checks on byte-identical inputs",
    "D06 accepted via exact equality with XVI-P P15 Sigma (ghost equation), not via the reviewer's own Feynman rules"]


def jl(p):
    return json.loads(Path(p).read_text())


def guard_ring():
    sys.path.insert(0, str(Q / "code" / "q01"))
    import qguard as G
    import sympy as sp
    from sympy import QQ
    from sympy.polys.rings import ring, sring
    from sympy.polys.fields import field
    d, X = G.d, sp.Symbol("X")
    R, D_, x = ring([d, X], QQ)
    F, Df, xf = field([d, X], QQ)
    forms = {"ring_evaluate": lambda: (D_**2 - D_ + x).evaluate(D_, 4), "ring_subs": lambda: (D_**2 - D_).subs(D_, 4),
             "ring_call": lambda: (D_**2 - D_)(4), "ring_compose": lambda: (D_**2).compose(D_, R(4)),
             "field_evaluate": lambda: ((Df - 1) / (Df - 4 + xf)).evaluate(Df, 5),
             "field_subs": lambda: ((Df - 1) / Df).subs(Df, 4), "field_call": lambda: ((Df - 1) / Df)(4),
             "sring_evaluate": lambda: (lambda Rp: Rp[1].evaluate(Rp[0].gens[0], 4))(sring(d**2 - d, d))}
    res = {}
    for k, f in forms.items():
        try:
            f()
            res[k] = "BYPASS"
        except G.PrematureD4Substitution:
            res[k] = "raises"
    R2, x2 = ring([X], QQ)
    ok = {"d_free_ring_evaluate": str((x2**2).evaluate(x2, 3)), "cancel": str(sp.cancel((d**2 - 16) / (d - 4))),
          "factor": str(sp.factor(d**2 - d)), "together": str(sp.together(1 / d + 1 / X))}
    return {"forms": res, "n_forms": len(res), "missed": [k for k, v in res.items() if v != "raises"],
            "legitimate_algebra": ok, "entry_points": G.ENTRY_POINTS,
            "residual_disclosed": "raw dense-list evaluation (densetools dup_eval/dmp_eval_in on lists carrying no symbol "
                                  "information) cannot be attributed to d and is not guarded; no package code uses it"}


def record():
    o1 = jl(CERT / "Q02_O1.json")
    bi = jl(CERT / "Q02_bulk_identity.json")
    ct = bi["bulk_counterterms_explicit"]
    ev = ct["evanescent_p22"]
    rep = jl(IND)
    v = rep["O1_verdict"]
    rep_sha = hashlib.sha256(IND.read_bytes()).hexdigest()
    binding_ok = rep_sha == BOUND["report_sha256"] and v["decision"] == BOUND["decision"] and \
        rep["source"]["zip_sha256"] == BOUND["reviewed_zip_sha256"] and rep["source"]["tree_sha256"] == BOUND["reviewed_tree_sha256"]
    ev_ops = ["K'.c E2 (evanescent)", "cbar.c SE (evanescent)"] + [f"h.h evanescent: {k}" for k in sorted(ev["hh_sum"])]
    domain = {"p14r_4D_directions": sorted(ct["coefficients"]), "n_p14r_4D": ct["n_directions"],
              "p22_evanescent_operators": ev_ops, "n_p22_evanescent": len(ev_ops),
              "n_enlarged_domain": ct["n_directions"] + len(ev_ops),
              "p22_evanescent_coefficients": {"Kc_E2": ev["Kc_E2"], "ccb_SE": ev["ccb_SE"], "hh_sum": ev["hh_sum"]},
              "statement": "The counterterm domain used for O.1 item 5 is the ENLARGED domain: the 11 p14r 4D directions "
                           "(over Q(i)) PLUS the p22 evanescent operators (disclosed p22 contract revision, outside p14r). "
                           "With the declared QdS external momenta the local divergent part contains evanescent poles that "
                           "are cancelled only by the p22 evanescent operators; 'ST4a in the declared p14r domain' alone "
                           "holds only for 4D externals. sigma_0 closure including the evanescent operators: verified "
                           "(p20/p22 independent reviews)."}
    items = []
    for i in o1["items"]:
        j = {k: i[k] for k in ("item", "owner", "XVI_O_status", "evidence", "earned")}
        if i["item"] == "Demonstrate identity remainder cancellation / CT solution":
            j["evidence"] = ("ST4a: MS-bar CT in the ENLARGED domain = 11 p14r 4D directions over Q(i) (phi.phi = 0 earned) "
                             "PLUS the p22 evanescent operators (K'.c E2, cbar.c SE, 8 h.h evanescent structures), which alone "
                             "cancel the evanescent poles with d-dimensional (QdS) externals; ST4b: c.h defect ZERO exactly in d "
                             "(P19, x = 0), c.b = ghost equation (Sigma == D06), c.cbar* unreachable, c.phi parity; full "
                             "sigma_0 closure of the complete MS-bar CT incl. the evanescent operators in every gh+1 sector "
                             "(p20/p22); cohomology of the domain 0; no free-EOM equivalence used")
            j["evidence_corrected_in"] = "Q03 (Q02 review finding ST4a_domain_wording, medium)"
        items.append(j)
    reasons = {}
    for k, val in o1["closure_reasons_answered"].items():
        reasons[k] = {"status": "ANSWERED", "answer": val.removeprefix("ANSWERED: ")}
    r7, r8 = list(o1["closure_reasons_answered"])[6], list(o1["closure_reasons_answered"])[7]
    reasons[r7]["answer"] += " Local matching uses the ENLARGED domain (p14r 4D + p22 evanescent)."
    reasons[r8] = {"status": "NOT_ANSWERED_IN_XVI_O_TERMS (superseded prescription)",
                   "answer": "XVI-O's o70-sync D11 matter-loop kernel counterterm prescription is NOT the one used and is "
                             "not completed: both XVI-O kernel continuations fail on evanescent external legs (p22 finding). "
                             "O.1 is closed under XVI-P's declared DRED prescription (QdS loop momenta, Q4S vertex indices, "
                             "MS-bar, Feynman +i0, p22 evanescent counterterm contract revision). The bulk sector (hh, phi phi, "
                             "cbar c, K'.c, one-point) is evaluated, not only the D11 kernel; this does not complete XVI-O's "
                             "kernel prescription."}
    allx = all(i["earned"] for i in items)
    closed = allx and binding_ok
    return {"landmark": "XVI-Q", "release": "Q03",
            "O1_closure": "CLOSED (scoped)" if closed else "NOT_CLOSED", "O1_closure_confirmed": closed,
            "confirmation": {**BOUND, "report_file": "reproducibility/independent_rebuild_report_q02.json",
                             "report_sha256_now": rep_sha, "binding_ok": binding_ok,
                             "reviewer_conditions_for_the_label": v["conditions_for_the_label"],
                             "conditions_discharged_by": ["counterterm_domain (enlarged)", "prescription"]},
            "scope_verbatim_from_independent_review": v["scope"],
            "prescription": {"used": "XVI-P declared prescription: DRED (QdS loop momenta, Q4S vertex indices, opaque d via "
                                     "qguard), MS-bar, Feynman +i0, p22 evanescent counterterms; m = 0; frozen gauge fermion",
                             "not_used": "XVI-O o70-sync D11 kernel counterterm prescription (closure reason 8): superseded "
                                         "on evanescent legs, both kernel continuations fail there (p22); reason 8 is NOT "
                                         "marked answered"},
            "counterterm_domain": domain, "items": items, "all_earned": allx, "closure_reasons": reasons,
            "enlarge_rules": o1["enlarge_rules"], "downgrade_rules": o1["downgrade_rules"],
            "not_implied": NOT_IMPLIED, "caveats": CAVEATS,
            "XVI_O_ledger": "read-only; O.1 there stays NOT_CLOSED (unchanged)"}


def main():
    (CERT / "Q03_O1_record.json").write_text(json.dumps(record(), indent=1, sort_keys=False) + "\n")
    (CERT / "Q03_guard_ring.json").write_text(json.dumps(guard_ring(), indent=1) + "\n")
    print("written Q03_O1_record.json, Q03_guard_ring.json")


if __name__ == "__main__":
    main()
