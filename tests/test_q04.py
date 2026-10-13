"""XVI-Q Q04 tests: engine validation (live, fast), inventory, gate fail-closed on honesty / value edits."""
import json, shutil, subprocess, sys  # noqa: E401
from pathlib import Path
import pytest

Q = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Q / "code" / "q04"))


def jl(p):
    return json.loads((Q / p).read_text())


def test_inventory_and_vertex():
    import p6_inventory as P
    i = P.inventory()
    assert i["equals_XVI_P_P23_15"] and i["n_E3_tadpoles"] == 3 and i["power_counting"]["verified_on_all_enumerated_graphs"]
    assert P.vertex_hhhphiphi_check() and not P.vertex_hhhphiphi_check(mutate=True)


def test_engine_d11_live():
    sys.path.insert(0, str(Q / "code" / "q04"))
    import q04_certs as C
    r = C.task_d11()
    assert r["ratio_mine_over_XVI_P"] == "I" and r["constant_ratio"]


def test_status_honest():
    s = jl("certificates/Q04_P6_status.json")
    assert s["P6_outcome"] == "Partial" and s["stages"]["counterterm_solve_kappa3"] == "Missing"
    e = jl("certificates/Q04_P6_eval.json")
    assert e["n_missing"] == 13 and e["n_evaluated_pole"] == 2 and e["n_derived_zero"] == 3


def _copy(tmp):
    dst = tmp / "XVI-Q"
    top = {"build", "dist", "archival", "manuscript"}
    shutil.copytree(Q, dst, ignore=lambda d, names: [n for n in names if n == "__pycache__" or (Path(d) == Q and n in top)])
    return dst


def _edit(root, rel, fn):
    p = root / rel
    o = json.loads(p.read_text())
    fn(o)
    p.write_text(json.dumps(o, indent=1) + "\n")


MUT = {
    "claim_closure": lambda r: _edit(r, "certificates/Q04_P6_status.json", lambda o: o["stages"].update(sigma0_closure_kappa3="closed")),
    "missing_to_zero": lambda r: _edit(r, "certificates/Q04_P6_eval.json", lambda o: o["graph_table"][0].update(status="derived 0")),
    "ward_flip": lambda r: _edit(r, "certificates/Q04_P6_eval.json", lambda o: o["matter_ward_three_point"]["declared"].update(residue_d4_zero=False)),
    "d11_ratio": lambda r: _edit(r, "certificates/Q04_P6_eval.json", lambda o: o["validation_D11"].update(ratio_mine_over_XVI_P="2*I")),
    "outcome_full": lambda r: _edit(r, "certificates/Q04_P6_status.json", lambda o: o.update(P6_outcome="Constructed")),
    "code_edit": lambda r: (r / "code/q04/uvpole.py").write_text((r / "code/q04/uvpole.py").read_text() + "\n#x\n"),
}


@pytest.mark.parametrize("name", sorted(MUT))
def test_gate_rejects(tmp_path, name):
    root = _copy(tmp_path)
    MUT[name](root)
    r = subprocess.run([sys.executable, "code/q04/release_q04.py"], cwd=root, capture_output=True, text=True)
    assert r.returncode == 1, r.stdout
    assert not (root / "RELEASE_STATUS.json").exists()


def test_gate_accepts(tmp_path):
    r = subprocess.run([sys.executable, "code/q04/release_q04.py"], cwd=_copy(tmp_path), capture_output=True, text=True)
    assert r.returncode == 0, r.stdout


def test_guard_scan_and_pins():
    import guard_scan as g
    r = g.report()
    assert r["scan_clean"] and r["probes"]["python_term_walk"].startswith("BYPASS")
    exp = json.loads((Q / "code/q04/q04_expected.json").read_text())
    assert "tools/source_set.py" in exp["file_sha256"]


def test_gate_rejects_source_set_edit(tmp_path):
    root = _copy(tmp_path)
    (root / "tools/source_set.py").write_text((root / "tools/source_set.py").read_text() + "\n#x\n")
    r = subprocess.run([sys.executable, "code/q04/release_q04.py"], cwd=root, capture_output=True, text=True)
    assert r.returncode == 1


def test_gate_rejects_escape_primitive(tmp_path):
    root = _copy(tmp_path)
    (root / "code/q04/uvpole.py").write_text((root / "code/q04/uvpole.py").read_text() + "\n_x = __import__('sympy').srepr(1)\n")
    r = subprocess.run([sys.executable, "code/q04/release_q04.py"], cwd=root, capture_output=True, text=True)
    assert r.returncode == 1
