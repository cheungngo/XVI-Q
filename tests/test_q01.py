"""XVI-Q Q01 tests: intake pins, O.2 items (values), guard forms, gate value checks (fail closed, no stale status)."""
from __future__ import annotations
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

Q = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Q / "code" / "q01"))
import q01_conventions as QC  # noqa: E402

CERT = Q / "certificates"


def cj(n):
    return json.loads((CERT / f"{n}.json").read_text())


def test_intake_pins_and_errata():
    it = cj("Q01_intake")
    assert it["pins_ok"] and it["XVI_P_p23_zip_ok"] and it["XVI_P_p23_extracted_identical"] and it["O1_ledger_sha_ok"]
    assert [e["id"] for e in it["errata"]][:5] == ["E1", "E2", "E3", "E4", "E5"]
    assert "overclaim" in it["errata"][0]["correction"].lower()


def test_O2_items_constructed_with_values():
    o2 = cj("Q01_O2")
    assert o2["item_status"] == QC.EXPECTED["O2_items"] and o2["O2_kappa2_master_status"] == "Constructed"
    a = o2["items"]["a_ell4"]["detail"]
    assert a["n_terms_mine"] == a["n_terms_export"] == 3360 and a["full_polynomial_equal"] and a["normalisation_from_XVI_O_probe"] == "-2"
    assert o2["W_terms"] == QC.EXPECTED["W_terms"]
    assert "NOT quantum" in o2["scope"]


def test_O2b_live_independent_engine():
    import o2b_antibracket as AB
    r = AB.main(Q / "build" / "xvip_run" / "exports")
    assert all(r[k]["antibracket_af0_vanishes"] for k in ("kappa0", "kappa1", "kappa2"))
    assert r["kappa2"]["W_terms"] == 15428


def test_master_mutations_all_break():
    m = cj("Q01_O2")["master_mutations"]
    assert {k: v["n_residuals"] for k, v in m.items()} == QC.EXPECTED["master_mutation_residuals"]
    assert all(v["claim_status"] != "verified" for v in m.values())


FORMS = ["Matrix_Dummy", "Matrix_y", "Poly_shift", "apart_list", "apart_residue", "captured_rem", "compose", "evalf", "groebner", "lambdify", "limit",
         "polys_ns_rem", "reduced", "rem", "rep_eval", "resultant", "subs", "xreplace"]


def test_guard_probe_all_forms_raise():
    g = cj("Q01_guard")
    assert g["forms"] == sorted(FORMS) and g["missed"] == [] and g["n_exprs"] == 4
    assert g["algebra"]["authorised_value"] == "1/80" and g["algebra"]["domain_only_poly"] == "True"


def test_guard_live_p23_bypasses():
    code = ("import sys; sys.path.insert(0, %r)\n" % str(Q / "code" / "q01") +
            "from sympy.polys.polytools import rem as r0\nimport qguard, sympy as sp, sympy.polys as spp\n"
            "d = qguard.d; out = []\n"
            "for f in (lambda: spp.rem(d**2-1, d-4, d), lambda: r0(d**2-1, d-4, d), lambda: sp.reduced(d**2-1, [d-4], d),\n"
            "          lambda: sp.groebner([d-4, d**2-1], d), lambda: sp.Poly(d**2-1, d).rep.eval(4),\n"
            "          lambda: sp.Matrix([1/(d**2-1)]).subs(d, sp.Dummy()).subs(sp.Dummy(), 4)):\n"
            "    try:\n        f(); out.append('NOT')\n    except qguard.PrematureD4Substitution:\n        out.append('R')\n"
            "print(','.join(out))")
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=300)
    assert r.stdout.strip() == "R,R,R,R,R,R", r.stdout + r.stderr


def test_release_status_flags_and_hard_stops():
    s = json.loads((CERT / "Q01_release.json").read_text())      # Q01 layer record (RELEASE_STATUS is Q02's)
    assert s["publication"] == "accepted" and s["execution_pass"] is True and all(v is True for v in s["value_checks"].values())
    assert s["O1_closure"] == "NOT_CLOSED" and s["bulk_symmetry_restored"] is False
    assert s["flags"]["O2_kappa2_master_constructed"] == "earned" and s["flags"]["bulk_O5_evaluated"] == "not_earned"
    assert s["hard_stops"] == QC.HARD_STOPS


def _copy():
    d = Path(tempfile.mkdtemp(prefix="q01gate_"))
    for sub in ("code", "certificates"):
        shutil.copytree(Q / sub, d / sub, ignore=shutil.ignore_patterns("__pycache__"))
    for f in ("RELEASE_STATUS.json",):
        if (Q / f).exists():
            shutil.copy(Q / f, d / f)
    return d


def _edit(d, n, fn):
    f = d / "certificates" / f"{n}.json"
    c = json.loads(f.read_text())
    fn(c)
    f.write_text(json.dumps(c))


MUTS = {
    "ell4_terms": lambda d: _edit(d, "Q01_O2", lambda c: c["items"]["a_ell4"]["detail"].__setitem__("n_terms_mine", 3359)),
    "ell4_norm": lambda d: _edit(d, "Q01_O2", lambda c: c["items"]["a_ell4"]["detail"].__setitem__("normalisation_from_XVI_O_probe", "2")),
    "W2_terms": lambda d: _edit(d, "Q01_O2", lambda c: c["W_terms"].__setitem__("kappa2", 15427)),
    "Ec_nonzero": lambda d: _edit(d, "Q01_O2", lambda c: c["items"]["b_antibracket_expansion"]["detail"]["unmutated"]["kappa2"]["E_c_W_terms"].__setitem__("0", 3)),
    "mut_count": lambda d: _edit(d, "Q01_O2", lambda c: c["master_mutations"]["corrupt_a2"].__setitem__("n_residuals", 0)),
    "status_c": lambda d: _edit(d, "Q01_O2", lambda c: c["item_status"].__setitem__("c_ghost_antifield", "Partial")),
    "export_sha": lambda d: _edit(d, "Q01_O2", lambda c: c["exports_sha256"].__setitem__("4", "0" * 64)),
    "guard_missed": lambda d: _edit(d, "Q01_guard", lambda c: c.__setitem__("missed", [["rem", "pow", "NOT_RAISED"]])),
    "pins": lambda d: _edit(d, "Q01_intake", lambda c: c.__setitem__("pins_ok", False)),
    "errata_dropped": lambda d: _edit(d, "Q01_intake", lambda c: c["errata"].pop(0)),
    "cert_missing": lambda d: (d / "certificates" / "Q01_O2.json").unlink(),
}


@pytest.mark.parametrize("name", sorted(MUTS))
def test_gate_rejects_value_edit_without_stale_status(name):
    d = _copy()
    try:
        MUTS[name](d)
        r = subprocess.run([sys.executable, str(d / "code" / "q01" / "release_q01.py")], capture_output=True, text=True)
        assert r.returncode == 1, r.stdout
        assert not (d / "certificates" / "Q01_release.json").exists()
    finally:
        shutil.rmtree(d, ignore_errors=True)


def test_gate_accepts_unmodified_copy():
    d = _copy()
    try:
        r = subprocess.run([sys.executable, str(d / "code" / "q01" / "release_q01.py")], capture_output=True, text=True)
        assert r.returncode == 0, r.stdout
    finally:
        shutil.rmtree(d, ignore_errors=True)
