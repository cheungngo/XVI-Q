#!/usr/bin/env python3
"""XVI-Q Q11 tests: gate fail-closed on value / honesty edits, unpinned files, and Q09-gate edits."""
import json, shutil, subprocess, sys  # noqa: E401
from pathlib import Path
import pytest

Q = Path(__file__).resolve().parents[7]
GATE = "code/p6q06/q07/q08/q09/q10/q11/release_q11.py"


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


SI, CL, S = ("certificates/Q11_sigma_images.json", "certificates/Q11_classification.json",
              "certificates/Q11_P6_status.json")
MUT = {
    "F_inflated": lambda r: _edit(r, SI, lambda o: o["results_at_derived_signs"]["F"]["derived_minus1"].update(in_N=194)),
    "reading_softened": lambda r: _edit(r, SI, lambda o: o.update(reading="all sigma-images close")),
    "unclassified_hidden": lambda r: _edit(r, CL, lambda o: o.update(unclassified=0)),
    "unconditional": lambda r: _edit(r, CL, lambda o: o.update(status="ESTABLISHED")),
    "physical_claimed": lambda r: _edit(r, CL, lambda o: o.update(physical_evanescent_directions="43 physical")),
    "basis_claimed": lambda r: _edit(r, S, lambda o: o["goals"].update(**{"3_invariant_basis": "DONE"})),
    "p22_claimed": lambda r: _edit(r, S, lambda o: o["goals"].update(**{"3_p22_consistency": "consistent"})),
    "old_label": lambda r: _edit(r, S, lambda o: o.update(sigma0_status_kappa3="CLOSED at pole level")),
    "insertion_unfrozen": lambda r: _edit(r, S, lambda o: o.update(insertion="open")),
    "revision_hidden": lambda r: _edit(r, S, lambda o: o.update(disclosed_revisions={})),
    "q04_report_edit": lambda r: (r / "reproducibility/independent_rebuild_report_q04.json").write_text("{}\n"),
    "unpinned_file": lambda r: (r / "code/p6q06/q07/q08/q09/q10/q11/sneak.py").write_text("x = 1\n"),
    "q10_gate_edit": lambda r: (r / "code/p6q06/q07/q08/q09/q10/release_q10.py").write_text(
        (r / "code/p6q06/q07/q08/q09/q10/release_q10.py").read_text() + "\n"),
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
