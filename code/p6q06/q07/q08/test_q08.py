#!/usr/bin/env python3
"""XVI-Q Q08 tests: gate fail-closed on value / honesty edits and unpinned files; subst_vec / ansatz sanity."""
import json, shutil, subprocess, sys  # noqa: E401
from pathlib import Path
import pytest

Q = Path(__file__).resolve().parents[4]
GATE = "code/p6q06/q07/q08/release_q08.py"


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


G_, F_, S = "certificates/Q08_ghost_equation.json", "certificates/Q08_free_directions.json", "certificates/Q08_P6_status.json"
MUT = {
    "rho": lambda r: _edit(r, G_, lambda o: o["loop_poles"].update(rho_total="2")),
    "ghost_unsolvable_hidden": lambda r: _edit(r, G_, lambda o: o.update(dim_free_with_ghost_eq=152)),
    "particular_unverified": lambda r: _edit(r, G_, lambda o: o["particular_solution"].update(verified_A_x_eq_b=False)),
    "all_trivial_claim": lambda r: _edit(r, F_, lambda o: o.update(dim_N_cap_all_families_upper_bound=152)),
    "fitted_hidden": lambda r: _edit(r, F_, lambda o: o.update(reading="all trivial")),
    "basis_claimed": lambda r: _edit(r, S, lambda o: o["goals"].update({"2_invariant_operator_basis": "DONE"})),
    "p22_claimed": lambda r: _edit(r, S, lambda o: o["goals"].update({"4_p22_consistency": "consistent"})),
    "revision_hidden": lambda r: _edit(r, S, lambda o: o.update(disclosed_revisions={})),
    "q07_value": lambda r: _edit(r, "certificates/Q07_hhc_row.json", lambda o: o.update(residual_4D="x")),
    "unpinned_code_file": lambda r: (r / "code/p6q06/q07/q08/extra.py").write_text("x = 1\n"),
    "q07_gate_edit": lambda r: (r / "code/p6q06/q07/release_q07.py").write_text((r / "code/p6q06/q07/release_q07.py").read_text() + "\n"),
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


def test_tmap_trace():
    sys.path.insert(0, str(Q / "code" / "p6q06" / "q07" / "q08"))
    import geq
    import uvpole as U
    assert geq._trace(U.H("u", "w")) == 4 and geq._trace(U.T("u", "w")) == 0
    assert geq._trace(U.H("u", "x") * U.H("w", "y")) == U.H("x", "y")
