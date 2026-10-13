#!/usr/bin/env python3
"""XVI-Q Q07 tests: gate fail-closed on value / honesty edits and unpinned files; subst_vec / ansatz sanity."""
import json, shutil, subprocess, sys  # noqa: E401
from pathlib import Path
import pytest

Q = Path(__file__).resolve().parents[3]
GATE = "code/p6q06/q07/release_q07.py"


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


def _gate(root):
    return subprocess.run([sys.executable, GATE], cwd=root, capture_output=True, text=True)


H, E, S = "certificates/Q07_hhc_row.json", "certificates/Q07_evanescent_CT.json", "certificates/Q07_P6_status.json"
MUT = {
    "hhc_4D_residual": lambda r: _edit(r, H, lambda o: o.update(residual_4D="x")),
    "sigma_F": lambda r: _edit(r, H, lambda o: o["tree_calibration"].update(sigma_F="-I")),
    "hhc_mutation_undetected": lambda r: _edit(r, H, lambda o: o["mutations_detected"].update(drop_E=False)),
    "remainder_value": lambda r: _edit(r, H, lambda o: o.update(evanescent_remainder="0")),
    "naive_claimed": lambda r: _edit(r, E, lambda o: o["naive_p22_type"]["ccK"].update(zero=True)),
    "CT_value": lambda r: _edit(r, E, lambda o: o["joint_local_CT"].update(particular_CT_Gamma_hhh="0")),
    "CT_unique_claim": lambda r: _edit(r, E, lambda o: o.update(verdict="CANCELLABLE uniquely")),
    "closure_overclaim": lambda r: _edit(r, S, lambda o: o.update(sigma0_closure_kappa3="CLOSED")),
    "not_implied_dropped": lambda r: _edit(r, S, lambda o: o.update(not_implied=[])),
    "revision_hidden": lambda r: _edit(r, S, lambda o: o.update(disclosed_revisions={})),
    "unpinned_code_file": lambda r: (r / "code/p6q06/q07/extra.py").write_text("x = 1\n"),
    "q06_gate_edit": lambda r: (r / "code/p6q06/release_q06.py").write_text((r / "code/p6q06/release_q06.py").read_text() + "\n"),
}


@pytest.mark.parametrize("name", sorted(MUT))
def test_gate_rejects(tmp_path, name):
    root = _copy(tmp_path)
    MUT[name](root)
    r = _gate(root)
    assert r.returncode != 0
    assert not (root / "RELEASE_STATUS.json").exists()


def test_gate_accepts(tmp_path):
    r = _gate(_copy(tmp_path))
    assert r.returncode == 0, r.stdout


def test_subst_vec_and_ansatz():
    sys.path.insert(0, str(Q / "code" / "p6q06" / "q07"))
    import hhc
    import evfit
    import uvpole as U
    e = U.H("gb2", "p1") + U.T("gb2", "a2")
    r = hhc.subst_vec(e, "gb2", {"a3": 2, "p1": 1})
    assert r == 2 * U.H("a3", "p1") + U.H("p1", "p1") + 2 * U.T("a3", "a2")
    ms = evfit.monomials({"z": 1, "gc1": 1}, 2)
    assert ms and all(any(s.name.startswith("T_") for s in m.free_symbols) for m in ms)
