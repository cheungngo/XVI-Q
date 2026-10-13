"""XVI-Q Q03 tests: O.1 record content, ring/field guard (live), gate fail-closed on corrections/provenance edits."""
import json, shutil, subprocess, sys  # noqa: E401
from pathlib import Path
import pytest

Q = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Q / "code" / "q01"))


def jl(p):
    return json.loads((Q / p).read_text())


def test_record():
    r = jl("certificates/Q03_O1_record.json")
    assert r["O1_closure"] == "CLOSED (scoped)" and r["O1_closure_confirmed"] is True and r["confirmation"]["binding_ok"]
    assert r["counterterm_domain"]["n_enlarged_domain"] == 21
    assert not list(r["closure_reasons"].values())[7]["status"].startswith("ANSWERED")
    assert len(r["caveats"]) == 4 and len(r["not_implied"]) == 6


def test_ring_field_guard_live():
    import qguard as G
    import sympy as sp
    from sympy import QQ
    from sympy.polys.rings import ring
    from sympy.polys.fields import field
    R, D_, x = ring([G.d, sp.Symbol("X")], QQ)
    F, Df, _ = field([G.d, sp.Symbol("X")], QQ)
    for f in (lambda: (D_ ** 2).evaluate(D_, 4), lambda: (D_ - 1)(4), lambda: ((Df - 1) / Df).subs(Df, 4)):
        with pytest.raises(G.PrematureD4Substitution):
            f()
    assert sp.cancel((G.d ** 2 - 16) / (G.d - 4)) == G.d + 4


def _copy(tmp):
    dst = tmp / "XVI-Q"
    top = {"build", "dist", "archival", "manuscript"}
    shutil.copytree(Q, dst, ignore=lambda d, names: [n for n in names if n == "__pycache__" or (Path(d) == Q and n in top)])
    (dst / "build" / "xvip_run").mkdir(parents=True)
    return dst


def _gate(root):
    return subprocess.run([sys.executable, "code/q03/release_q03.py"], cwd=root, capture_output=True, text=True)


def _edit_json(root, rel, fn):
    p = root / rel
    o = json.loads(p.read_text())
    fn(o)
    p.write_text(json.dumps(o, indent=1) + "\n")


MUT = {
    "provenance_P19_edit": lambda r: (r / "provenance/XVI-P/certificates/P19_defect.json").write_text(
        (r / "provenance/XVI-P/certificates/P19_defect.json").read_text().replace("ZERO", "NONZERO", 1)),
    "unpinned_provenance_file": lambda r: (r / "provenance/XVI-O/extra.json").write_text("{}"),
    "reason8_marked_answered": lambda r: _edit_json(r, "certificates/Q03_O1_record.json", lambda o: list(
        o["closure_reasons"].values())[7].update(status="ANSWERED")),
    "domain_shrunk": lambda r: _edit_json(r, "certificates/Q03_O1_record.json", lambda o: o["counterterm_domain"].update(n_p22_evanescent=0)),
    "caveat_dropped": lambda r: _edit_json(r, "certificates/Q03_O1_record.json", lambda o: o["caveats"].pop(1)),
    "report_edited": lambda r: (r / "reproducibility/independent_rebuild_report_q02.json").write_text(
        (r / "reproducibility/independent_rebuild_report_q02.json").read_text() + " "),
    "qguard_edited": lambda r: (r / "code/q01/qguard.py").write_text((r / "code/q01/qguard.py").read_text() + "\n#x\n"),
    "expected_edited": lambda r: (r / "code/q03/q03_expected.json").write_text((r / "code/q03/q03_expected.json").read_text() + " "),
}


@pytest.mark.parametrize("name", sorted(MUT))
def test_gate_rejects(tmp_path, name):
    root = _copy(tmp_path)
    MUT[name](root)
    r = _gate(root)
    assert r.returncode == 1, r.stdout
    assert not (root / "RELEASE_STATUS.json").exists()


def test_gate_accepts_unmodified(tmp_path):
    r = _gate(_copy(tmp_path))
    assert r.returncode == 0, r.stdout
