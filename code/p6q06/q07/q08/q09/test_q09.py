#!/usr/bin/env python3
"""XVI-Q Q09 tests: gate fail-closed on value / honesty edits and unpinned files; subst_vec / ansatz sanity."""
import json, shutil, subprocess, sys  # noqa: E401
from pathlib import Path
import pytest

Q = Path(__file__).resolve().parents[5]
GATE = "code/p6q06/q07/q08/q09/release_q09.py"


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


T_, S, B_ = "certificates/Q09_trivial_families.json", "certificates/Q09_P6_status.json", "certificates/Q09_q04_report_binding.json"
MUT = {
    "lambda_F": lambda r: _edit(r, T_, lambda o: o["families"]["F"].update(lambda_derived="1")),
    "F_all_in_N": lambda r: _edit(r, T_, lambda o: o["families"]["F"].update(n_in_N=204, n_fail=0)),
    "span_inflated": lambda r: _edit(r, T_, lambda o: o.update(dim_trivial_span_derived=104)),
    "old_label": lambda r: _edit(r, S, lambda o: o.update(sigma0_status_kappa3="CLOSED at pole level")),
    "basis_claimed": lambda r: _edit(r, S, lambda o: o["goals"].update(c_invariant_operator_basis="DONE")),
    "p22_claimed": lambda r: _edit(r, S, lambda o: o["goals"].update(d_p22_consistency="consistent")),
    "revision_hidden": lambda r: _edit(r, S, lambda o: o.update(disclosed_revisions={})),
    "q04_report_edit": lambda r: (r / "reproducibility/independent_rebuild_report_q04.json").write_text("{}\n"),
    "scan3_rules_edit": lambda r: (r / "code/p6q06/q07/q08/q09/guard_scan3_rules.json").write_text('{"dict_walk_methods": []}\n'),
    "methodcaller_snuck_in": lambda r: (r / "code/p6q06/q07/q08/q09/sneak.py").write_text("from operator import methodcaller\n"),
    "q08_gate_edit": lambda r: (r / "code/p6q06/q07/q08/release_q08.py").write_text((r / "code/p6q06/q07/q08/release_q08.py").read_text() + "\n"),
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


def test_scan3_probes():
    sys.path.insert(0, str(Q / "code" / "p6q06" / "q07" / "q08" / "q09"))
    import ast
    import guard_scan3 as g
    assert all(g.hits(ast.parse(v)) for v in g.PROBES.values())
    assert not g.hits(ast.parse("e.coeff(t, 3)\nx.as_powers_dict()[s]"))
