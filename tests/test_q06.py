#!/usr/bin/env python3
"""XVI-Q Q06 tests: gate fail-closed on value / honesty edits and on unpinned code files; fast path sanity."""
import json, shutil, subprocess, sys  # noqa: E401
from pathlib import Path
import pytest

Q = Path(__file__).resolve().parents[1]


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
    return subprocess.run([sys.executable, "code/p6q06/release_q06.py"], cwd=root, capture_output=True, text=True)


def _gr(o, t, **kw):
    o[t].update(**kw)


MUT = {
    "row_4D_residual": lambda r: _edit(r, "certificates/Q06_ghost_rows.json", lambda o: _gr(o, "c.c.K", residual_4D="x")),
    "evanescent_claimed_closed": lambda r: _edit(r, "certificates/Q06_ghost_rows.json", lambda o: o.update(status="CLOSED")),
    "K_sign_back_to_Q05": lambda r: _edit(r, "certificates/Q06_ghost_rows.json",
                                          lambda o: o["c.c.cbar"]["channel_coefficients"].update(K1="1")),
    "mutation_undetected": lambda r: _edit(r, "certificates/Q06_ghost_rows.json",
                                           lambda o: o["c.c.K"]["mutations_detected"].update(drop_c1=False)),
    "V4_gauge": lambda r: _edit(r, "certificates/Q06_V4_validation.json", lambda o: o.update(gauge_n3_V4_vs_V3_zero=False)),
    "G01_value": lambda r: _edit(r, "certificates/Q06_P6_hhh.json",
                                 lambda o: o["graphs"]["G01"].update(pole_over_U_d4_residue="0")),
    "claim_closure": lambda r: _edit(r, "certificates/Q06_P6_status.json", lambda o: o.update(sigma0_closure_kappa3="closed")),
    "hhc_claimed": lambda r: _edit(r, "certificates/Q06_P6_status.json", lambda o: o["stages"].update(row_h_h_c="closed")),
    "unpinned_code_file": lambda r: (r / "code/p6q06/extra.py").write_text("x = 1\n"),
    "q05_gate_edit": lambda r: (r / "code/p6q05/release_q05.py").write_text((r / "code/p6q05/release_q05.py").read_text() + "\n"),
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


def test_fast_trunc_exact():
    sys.path.insert(0, str(Q / "code" / "p6q06"))
    import fastpole as FP
    import uvpole as U
    F = FP.Fast(["p1", "a1"])
    e = U.t ** 2 * U.H("a1", "p1") + U.t ** 5 * U.H("p1", "p1")
    assert F.back(F.trunc(F.to(e), 4)) == U.t ** 2 * U.H("a1", "p1")
