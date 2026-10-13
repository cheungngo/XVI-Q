#!/usr/bin/env python3
"""XVI-Q Q01: writes certificates Q01_intake.json, Q01_O2.json, Q01_guard.json from the build products
(build/xvip_run = re-executed XVI-P p23 classical stages; build/mut_* = master mutations)."""
from __future__ import annotations
import hashlib
import json
import subprocess
import sys
from pathlib import Path

Q = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Q / "code" / "q01"))
import q01_conventions as QC  # noqa: E402
import ell4_independent as L4  # noqa: E402
import o2b_antibracket as AB  # noqa: E402

RUN = Q / "build" / "xvip_run"
CERT = Q / "certificates"


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def content_sha(p):
    sys.path.insert(0, str(RUN / "code" / "classical"))
    import certio
    return certio.file_content_sha(p, ())


def intake():
    imp = json.loads((Q / "intake" / "imports.json").read_text())
    ok = {k: (Q / k).is_file() and sha(Q / k) == v["sha256"] for k, v in imp["files"].items()}
    z = __import__("zipfile").ZipFile(Q / "provenance/XVI-P/XVI-P-p23-source-only.zip")
    extracted_ok = all(sha(Q / "imports" / ("XVI-P-p23" + n[len("XVI-P"):])) == hashlib.sha256(z.read(n)).hexdigest()
                       for n in z.namelist() if not n.endswith("/"))
    return {"pins_ok": all(ok.values()), "n_pins": len(ok), "pins": ok,
            "XVI_P_p23_zip_ok": imp["XVI_P_p23_zip_sha256"] == QC.XVI_P_P23_ZIP_SHA256,
            "XVI_P_p23_extracted_identical": extracted_ok, "n_extracted": len([n for n in z.namelist() if not n.endswith("/")]),
            "O1_ledger_sha_ok": imp["files"]["provenance/XVI-O/ledgers/O_quantum_target.json"]["sha256"] == QC.O1_LEDGER_SHA256,
            "errata": json.loads((Q / "intake" / "errata_p23.json").read_text())["errata"],
            "registration": json.loads((Q / "registration" / "local_registration.json").read_text())["registered"]}


def o2():
    ex = {k: sha(RUN / "exports" / f"ell{k}_generating_polynomial.srepr") for k in (2, 3, 4)}
    a = L4.main(RUN / "exports" / "ell4_generating_polynomial.srepr")
    b = {"unmutated": AB.main(RUN / "exports")}
    for m in ("scale_ell4_term", "drop_s1"):
        import importlib
        importlib.reload(AB)
        r = AB.main(RUN / "exports", mutate=m)
        b[f"mutation_{m}_detected"] = not r["kappa2"]["antibracket_af0_vanishes"]
    importlib.reload(AB)
    master = json.loads((RUN / "certificates" / "P13_master_kappa2.json").read_text())
    brst = json.loads((RUN / "certificates" / "P13_brst_derivation.json").read_text())
    same = {c: content_sha(RUN / "certificates" / f"{c}.json") == content_sha(Q / "provenance" / "XVI-P" / "certificates" / f"{c}.json")
            for c in QC.EXPECTED["classical_certs_content_identical_to_p23"]}
    muts = {}
    for m in QC.EXPECTED["master_mutation_residuals"]:
        c = json.loads((Q / "build" / f"mut_{m}" / "certificates" / "P13_master_kappa2.json").read_text())
        muts[m] = {"claim_status": c["claim_status"], "n_residuals": len(c["required_master_residuals"])}
    tbl = master["table_by_kappa_and_antifield_number"]["kappa2"]
    items = {
        "a_ell4": {"status": "Constructed" if (a["full_polynomial_equal"] and a["export_on_XVI_O_probe_equal"]
                                                and ex[4] == QC.EXPORT_SHA256[4]) else "Partial",
                   "argument": "XVI-P p23 ell_4 export (3360 terms, re-executed byte-identical) equals an INDEPENDENT "
                               "Gamma-Gamma reconstruction term by term; the single normalisation (-2) is fixed by "
                               "XVI-O's conformal probe and equals XVI-O's convention ell_4 = -2[eps^4](sqrt(-q)R); the "
                               "export also reproduces XVI-O's probe exactly", "detail": a},
        "b_antibracket_expansion": {
            "status": "Constructed" if all(b["unmutated"][k]["antibracket_af0_vanishes"] for k in ("kappa0", "kappa1", "kappa2"))
            and b["mutation_scale_ell4_term_detected"] and b["mutation_drop_s1_detected"] else "Partial",
            "argument": "(S,S)|af0,grav at kappa^2 = 2 int W_2, W_2 = sum E(ell_4) s0 h + E(ell_3) s1 h (s2 = 0 derived) -- the "
                        "explicit component form of 2(S0,S2)+(S1,S1); W_2 has 15428 terms and E_{c^mu}(W_2) = 0 exactly "
                        "(independent XVI-Q jet engine; agrees with XVI-P grav_noether re-executed)", "detail": b},
        "c_ghost_antifield": {
            "status": "Constructed" if (all(tbl.values()) and brst["claim_status"] == "s2_derived_zero_all_sectors"
                                        and brst["nilpotency_exact_all_orders"] and all(same.values())
                                        and all(muts[m]["n_residuals"] == n for m, n in QC.EXPECTED["master_mutation_residuals"].items()))
            else "Partial",
            "argument": "re-executed exact XVI-P P13 master: af1/af2 (invariant and Psi-transformed variables) and the af0 "
                        "ghost/gauge-fixing cell vanish at kappa^0..kappa^2 with s2 DERIVED zero in all sectors (not "
                        "invented; answers XVI-O's 'Absent under inherited s'); s nilpotent exactly; 7 mutations each "
                        "produce residuals. Evidence = re-execution of XVI-P's engine (no second independent engine for "
                        "this sector in Q01).", "kappa2_table": tbl, "s2": brst["s2_per_sector"]}}
    st = {k: v["status"] for k, v in items.items()}
    o2_status = "Constructed" if all(v == "Constructed" for v in st.values()) else "Partial"
    return {"equation": "2(S0,S2) + (S1,S1) = 0 on S_inv (classical BV master at kappa^2)",
            "items": items, "item_status": st, "O2_kappa2_master_status": o2_status,
            "scope": "classical, S_inv, through O(kappa^2), m = 0; NOT quantum, NOT the prepared S_D, NOT the insertion",
            "exports_sha256": ex, "exports_identical_to_p23": all(ex[k] == QC.EXPORT_SHA256[k] for k in ex),
            "classical_certs_content_identical_to_p23": same, "master_mutations": muts,
            "W_terms": {k: b["unmutated"][k]["W_terms"] for k in ("kappa0", "kappa1", "kappa2")}}


def guard():
    p = subprocess.run([sys.executable, str(Q / "code" / "q01" / "guard_probe.py")], capture_output=True, text=True,
                       timeout=600)
    return json.loads(p.stdout.strip().splitlines()[-1])


def main():
    CERT.mkdir(exist_ok=True)
    out = {"Q01_intake": intake(), "Q01_O2": o2(), "Q01_guard": guard()}
    for k, v in out.items():
        (CERT / f"{k}.json").write_text(json.dumps({"landmark": "XVI-Q", "release": "Q01", **v}, indent=1, default=str) + "\n")
    print(json.dumps({"intake": out["Q01_intake"]["pins_ok"], "O2": out["Q01_O2"]["item_status"],
                      "guard_missed": out["Q01_guard"]["missed"]}))


if __name__ == "__main__":
    main()
