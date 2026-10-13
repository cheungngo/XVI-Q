"""XVI-Q Q02 tests: bulk graphs, D12 reconciliation, bulk identity, O.1 walk, guard (apart), gate value checks."""
from __future__ import annotations
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

Q = Path(__file__).resolve().parents[1]
CERT = Q / "certificates"


def cj(n):
    return json.loads((CERT / f"{n}.json").read_text())


def test_D04_zero_by_two_routes_and_factorisation():
    d4 = cj("Q02_bulk_graphs")["D04"]
    assert d4["value"] == "0" and d4["routes_agree"] and d4["route1_A"] == "0" and d4["route2_IBP_A"] == "0"
    assert d4["N_proportional_to_p2_times_scalar_denominator"] is True


def test_D04_hand_contraction_live():
    import sympy as sp
    p, q, pq = sp.symbols("p2 q2 pq")
    # T.P.T' = 2(p.p')(q.q') + 2(p.q')(q.p') - 2(p.q)(p'.q'), p' = -p, q' = -q
    val = 2 * (-p) * (-q) + 2 * (-pq) * (-pq) - 2 * pq * pq
    assert sp.expand(val - 2 * p * q) == 0


def test_D06_equals_XVI_P_and_poles():
    d6 = cj("Q02_bulk_graphs")["D06"]
    assert d6["equals_XVI_P_P15_exactly"] and d6["fit_exact"]
    assert {k: v["UV_pole"] for k, v in d6["laurent"].items()} == {"S1_delta^rho_nu_p4": "-3/8", "S2_p^rho_p_nu_p2": "1/4",
                                                                   "SE_Gt^rho_nu_p4": "-1/8"}


def test_tadpoles_derived_no_log():
    tp = cj("Q02_bulk_graphs")["tadpoles"]
    assert tp["D05a"]["omegas"] == [2] and tp["D12"]["omegas"] == [2, 4]
    assert not tp["D05a"]["logarithmic_member"] and not tp["D12"]["logarithmic_member"]


def test_D12_reconciliation():
    r = cj("Q02_D12_reconciliation")
    assert r["consistent"] and r["XVI_P_p23_stale_line_present"] and "Missing" not in r["status"]


def test_bulk_identity_rows_and_closure():
    b = cj("Q02_bulk_identity")
    assert b["n_rows"] == 23 and b["unclassified_rows"] == [] and b["bulk_identity_outcome"] == "cancellable"
    assert b["rows"]["phi.phi"]["outcome"] == "zero" and b["rows"]["b.c"]["class"] == "ST4b_ghost_equation"
    assert b["sigma0_closure"]["closed"] and b["ST4a"]["cbc"]["SE_pole_plus_p22_CT_zero"]


def test_O1_walk_and_status():
    o = cj("Q02_O1")
    assert o["all_earned"] and o["O1_closure"] == "CLOSED_PENDING_INDEPENDENT_REVIEW" and len(o["items"]) == 6
    assert o["O1_closure_confirmed"] is False and len(o["closure_reasons_answered"]) == 8
    s = json.loads((Q / "certificates" / "Q02_release.json").read_text())   # Q02 layer record (RELEASE_STATUS is Q03's)
    assert s["O1_closure"] == "CLOSED_PENDING_INDEPENDENT_REVIEW" and s["O1_closure_confirmed"] is False
    assert s["bulk_symmetry_restored"] is False and "QdS" in s["O1_closure_scope"]


def test_D04_route4_live():
    sys.path.insert(0, str(Q / "code" / "q02"))
    import d04_route4 as R4
    r = R4.main()
    assert r["declared_zero"] and r["all_mutations_nonzero"] and r["T3_zero_by_Ward_identity"]
    assert r == {k: v for k, v in cj("Q02_D04_scrutiny").items() if k in r}


def test_enumeration_and_CTs():
    b = cj("Q02_bulk_identity")
    assert b["independent_graph_enumeration"]["n_E2"] == 9 and b["independent_graph_enumeration"]["n_E1"] == 3
    c = b["bulk_counterterms_explicit"]
    assert c["ok"] and c["coefficients"]["h.h#0"] == "49*I/480" and c["coefficients"]["phi.phi#0"] == "0"


def test_guard_apart_bypass_live():
    code = ("import sys; sys.path.insert(0, %r)\n" % str(Q / "code" / "q01") +
            "import qguard, sympy as sp\nd = qguard.d\n"
            "try:\n    sp.apart((d**2 - 1)/(d - 4), d); print('NOT')\nexcept qguard.PrematureD4Substitution:\n    print('R')\n")
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=300)
    assert r.stdout.strip() == "R", r.stdout + r.stderr


def _copy():
    d = Path(tempfile.mkdtemp(prefix="q02gate_"))
    for sub in ("code", "certificates", "intake", "registration"):
        shutil.copytree(Q / sub, d / sub, ignore=shutil.ignore_patterns("__pycache__"))
    (d / "build" / "xvip_run" / "certificates").mkdir(parents=True)
    shutil.copy(Q / "build/xvip_run/certificates/P13_master_kappa2.json", d / "build/xvip_run/certificates/")
    return d


def _edit(d, n, fn):
    f = d / "certificates" / f"{n}.json"
    c = json.loads(f.read_text())
    fn(c)
    f.write_text(json.dumps(c))


def _jedit(path, fn):
    c = json.loads(path.read_text())
    fn(c)
    path.write_text(json.dumps(c))


MUTS = {
    "O2c_cell_false": lambda d: _edit(d, "Q01_O2", lambda c: c["items"]["c_ghost_antifield"]["kappa2_table"].__setitem__("af2_transformed", False)),
    "O2c_table_replaced": lambda d: _edit(d, "Q01_O2", lambda c: c["items"]["c_ghost_antifield"].__setitem__("kappa2_table", {"x": True})),
    "O2c_build_input_cell": lambda d: _jedit(d / "build/xvip_run/certificates/P13_master_kappa2.json",
                                             lambda c: c["table_by_kappa_and_antifield_number"]["kappa1"].__setitem__("af1_invariant", False)),
    "per_file_pin_false": lambda d: _edit(d, "Q01_intake", lambda c: c["pins"].__setitem__(next(iter(c["pins"])), False)),
    "expected_file_edit": lambda d: (d / "code/q02/q02_expected.json").write_text((d / "code/q02/q02_expected.json").read_text().replace("-3/8", "-1/8", 1)),
    "D04_value": lambda d: _edit(d, "Q02_bulk_graphs", lambda c: c["D04"].__setitem__("route1_A", "1/2")),
    "D06_pole": lambda d: _edit(d, "Q02_bulk_graphs", lambda c: c["D06"]["laurent"]["S2_p^rho_p_nu_p2"].__setitem__("UV_pole", "1/3")),
    "D06_local_value": lambda d: _edit(d, "Q02_bulk_graphs", lambda c: c["D06"]["laurent"]["S1_delta^rho_nu_p4"].__setitem__("local", "-1/3")),
    "row_outcome": lambda d: _edit(d, "Q02_bulk_identity", lambda c: c["rows"]["phi.phi"].__setitem__("outcome", "obstruction")),
    "O1_item_unearned": lambda d: _edit(d, "Q02_O1", lambda c: c["items"][3].__setitem__("earned", False)),
    "errata_E2_severity": lambda d: _jedit(d / "intake/errata_p23.json", lambda c: c["errata"][1].__setitem__("severity", "info")),
    "registration_reverted": lambda d: _jedit(d / "registration/local_registration.json", lambda c: c.__setitem__("how_found", "old")),
    "hh_per_graph_ct_edit": lambda d: _edit(d, "Q02_bulk_identity", lambda c: c["bulk_counterterms"]["cbc"].__setitem__("SE", "1/4")),
    "D04_scrutiny_mutation_zeroed": lambda d: _edit(d, "Q02_D04_scrutiny", lambda c: c["cases"]["MUT non-de Donder gauge T4 piece (beta != -2), g4=1"].__setitem__("value_over_master", "0")),
    "O1_label_confirmed_early": lambda d: _edit(d, "Q02_O1", lambda c: c.__setitem__("O1_closure_confirmed", True)),
    "CT_phiphi_edit": lambda d: _edit(d, "Q02_bulk_identity", lambda c: c["bulk_counterterms_explicit"]["coefficients"].__setitem__("phi.phi#0", "1/2")),
    "conventions_edit": lambda d: (d / "code/q01/q01_conventions.py").write_text((d / "code/q01/q01_conventions.py").read_text() + "\n# edit\n"),
}


@pytest.mark.parametrize("name", sorted(MUTS))
def test_gate_rejects_without_stale_status(name):
    d = _copy()
    try:
        MUTS[name](d)
        r = subprocess.run([sys.executable, str(d / "code/q02/release_q02.py")], capture_output=True, text=True)
        assert r.returncode == 1, r.stdout
        assert not (d / "RELEASE_STATUS.json").exists() and not (d / "certificates/Q02_release.json").exists()
    finally:
        shutil.rmtree(d, ignore_errors=True)


def test_gate_accepts_unmodified_copy():
    d = _copy()
    try:
        r = subprocess.run([sys.executable, str(d / "code/q02/release_q02.py")], capture_output=True, text=True)
        assert r.returncode == 0, r.stdout
    finally:
        shutil.rmtree(d, ignore_errors=True)
