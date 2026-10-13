#!/usr/bin/env python3
"""XVI-Q Q05 tests: hardened scan, gate fail-closed on value / honesty edits and on unpinned code files."""
import json, shutil, subprocess, sys  # noqa: E401
from pathlib import Path
import pytest

Q = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Q / "code" / "p6q05"))


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
    return subprocess.run([sys.executable, "code/p6q05/release_q05.py"], cwd=root, capture_output=True, text=True)


MUT = {
    "row_residual": lambda r: _edit(r, "certificates/Q05_P6_rows.json", lambda o: o["c.phi.phi"].update(residual="x")),
    "ghost_row_closed": lambda r: _edit(r, "certificates/Q05_P6_rows.json",
                                        lambda o: o["ghost_rows"]["c.c.K"].update(loop_residual_zero=True)),
    "claim_closure": lambda r: _edit(r, "certificates/Q05_P6_status.json",
                                     lambda o: o.update(sigma0_closure_kappa3="closed")),
    "missing_to_zero": lambda r: _edit(r, "certificates/Q05_P6_graphs.json", lambda o: o["missing"].pop("G01")),
    "D13_ratio": lambda r: _edit(r, "certificates/Q05_rules_validation.json",
                                 lambda o: o["validations"]["D13"].update(ratio="I")),
    "pole_value": lambda r: _edit(r, "certificates/Q05_P6_graphs.json",
                                  lambda o: o["graphs"]["G12"].update(pole_over_U_d4_residue="0")),
    "unpinned_code_file": lambda r: (r / "code/p6q05/sub").mkdir() or (r / "code/p6q05/sub/extra.py").write_text("x = 1\n"),
    "unpinned_tool": lambda r: (r / "tools/extra_tool.py").write_text("x = 1\n"),
    "escape_primitive": lambda r: (r / "code/p6q05/graphs.py").write_text((r / "code/p6q05/graphs.py").read_text()
                                                                         + "\nfrom sympy import lambdify as _f\n"),
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


def test_scan2_flags_reviewer_forms():
    import guard_scan2 as g
    r = g.report()
    assert r["all_reviewer_evasions_flagged"] and r["scan_clean"] and "NOT a security boundary" in r["claim"]
