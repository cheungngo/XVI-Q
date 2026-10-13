#!/usr/bin/env python3
"""XVI-Q Q10 tests: gate fail-closed on value / honesty edits, unpinned files, and Q09-gate edits."""
import json, shutil, subprocess, sys  # noqa: E401
from pathlib import Path
import pytest

Q = Path(__file__).resolve().parents[6]
GATE = "code/p6q06/q07/q08/q09/q10/release_q10.py"


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


A_, E_, S = ("certificates/Q10_family_anatomy.json", "certificates/Q10_erratum_Q08_counts.json",
             "certificates/Q10_P6_status.json")
MUT = {
    "blind_count": lambda r: _edit(r, A_, lambda o: o["blind_coordinates"].update(n=0)),
    "F_lambda_nonvacuous": lambda r: _edit(r, A_, lambda o: o["families"]["F"].update(lam_any=0)),
    "G_fixed": lambda r: _edit(r, A_, lambda o: o["families"]["G"].update(lam_none=0)),
    "erratum_dropped": lambda r: _edit(r, E_, lambda o: o["withdrawn"].update(Q08=[])),
    "caveat_dropped": lambda r: _edit(r, E_, lambda o: o["restated_bound"].update(caveat="")),
    "basis_claimed": lambda r: _edit(r, S, lambda o: o["goals"].update(invariant_evanescent_basis="DONE")),
    "GP_claimed": lambda r: _edit(r, S, lambda o: o["goals"].update(G_P_families="DONE")),
    "old_label": lambda r: _edit(r, S, lambda o: o.update(sigma0_status_kappa3="CLOSED at pole level")),
    "revision_hidden": lambda r: _edit(r, S, lambda o: o.update(disclosed_revisions={})),
    "intake_edit": lambda r: _edit(r, "intake/new_manuscripts_q10.json", lambda o: o.update(premises_imported=["En05-B"])),
    "q04_report_edit": lambda r: (r / "reproducibility/independent_rebuild_report_q04.json").write_text("{}\n"),
    "unpinned_file": lambda r: (r / "code/p6q06/q07/q08/q09/q10/sneak.py").write_text("x = 1\n"),
    "q09_gate_edit": lambda r: (r / "code/p6q06/q07/q08/q09/release_q09.py").write_text(
        (r / "code/p6q06/q07/q08/q09/release_q09.py").read_text() + "\n"),
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
