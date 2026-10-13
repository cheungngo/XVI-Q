#!/usr/bin/env python3
"""XVI-Q Q07 rebuild entrypoint.  --write-reference (author, built tree): reproducibility/reference_values_q07.json =
content sha256 (run-specific keys 'seconds' removed) of every Q01/Q02 certificate + RELEASE_STATUS + export shas.
Default: copy the source-only zip's tree to a temp dir, run tools/build.py there (no sync), compare with the reference."""
from __future__ import annotations
import hashlib, json, os, shutil, subprocess, sys, tempfile, zipfile  # noqa: E401
from pathlib import Path

Q = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Q / "tools"))
import source_set_q06_q07 as SS  # noqa: E402
REF = "reproducibility/reference_values_q07.json"
CERTS = ["certificates/Q01_intake.json", "certificates/Q01_O2.json", "certificates/Q01_guard.json",
         "certificates/Q01_release.json", "RELEASE_STATUS.json", "ledgers/Q01_milestones.json",
         "certificates/Q02_bulk_graphs.json", "certificates/Q02_D12_reconciliation.json", "certificates/Q02_bulk_identity.json",
         "certificates/Q02_O1.json", "certificates/Q02_D04_scrutiny.json", "certificates/Q02_release.json", "ledgers/Q02_milestones.json",
         "certificates/Q03_O1_record.json", "certificates/Q03_guard_ring.json", "certificates/Q03_release.json", "ledgers/Q03_milestones.json",
         "certificates/Q04_P6_inventory.json", "certificates/Q04_P6_identity.json", "certificates/Q04_P6_eval.json",
         "certificates/Q04_P6_status.json", "certificates/Q04_release.json", "ledgers/Q04_milestones.json",
         "certificates/Q04_guard_scan.json", "certificates/Q05_rules_validation.json", "certificates/Q05_P6_graphs.json",
         "certificates/Q05_P6_rows.json", "certificates/Q05_P6_status.json", "certificates/Q05_guard_scan.json",
         "certificates/Q05_release.json", "ledgers/Q05_milestones.json",
         "certificates/Q06_ghost_rows.json", "certificates/Q06_G10_crosscheck.json", "certificates/Q06_V4_validation.json",
         "certificates/Q06_fastpath.json", "certificates/Q06_P6_hhh.json", "certificates/Q06_P6_status.json",
         "certificates/Q06_release.json", "certificates/Q07_hhc_row.json", "certificates/Q07_evanescent_CT.json",
         "certificates/Q07_P6_status.json", "certificates/Q07_release.json", "ledgers/Q07_milestones.json"]


def strip(o):
    if isinstance(o, dict):
        return {k: strip(v) for k, v in o.items() if k != "seconds"}
    if isinstance(o, list):
        return [strip(v) for v in o]
    return o


def values(root):
    v = {c: hashlib.sha256(json.dumps(strip(json.loads((root / c).read_text())), sort_keys=True).encode()).hexdigest()
         for c in CERTS}
    v.update({f"export_ell{k}": hashlib.sha256((root / f"build/xvip_run/exports/ell{k}_generating_polynomial.srepr")
                                               .read_bytes()).hexdigest() for k in (2, 3, 4)})
    return v


def main():
    if "--write-reference" in sys.argv:
        (Q / "reproducibility").mkdir(exist_ok=True)
        (Q / REF).write_text(json.dumps({"tree_sha256_without_reference": SS.tree_sha(Q, exclude=(REF,)),
                                         "values": values(Q)}, indent=1) + "\n")
        print("written", REF)
        return 0
    ref = json.loads((Q / REF).read_text())
    with tempfile.TemporaryDirectory(prefix="q07src_") as td:
        zp = Path(td) / "src.zip"
        SS.make_zip(zp, Q)
        zipfile.ZipFile(zp).extractall(td)
        root = Path(td) / "XVI-Q"
        r = subprocess.run([sys.executable, "tools/build_q06_q07.py"], cwd=root, capture_output=True, text=True,
                           env=dict(os.environ, Q_SYNC_DIR="none"))
        (Q / "build").mkdir(exist_ok=True)
        (Q / "build" / "rebuild_build.log").write_text(r.stdout + r.stderr)
        built = r.returncode == 0 and "BUILD_OK" in r.stdout
        got = values(root) if built else {}
        diff = {k: (v, got.get(k)) for k, v in ref["values"].items() if got.get(k) != v}
        v_files = {f: hashlib.sha256((Q / f).read_bytes()).hexdigest() for f in ('code/q02/release_q02.py', 'code/q02/q02_expected.json', 'code/q03/release_q03.py', 'code/q03/q03_expected.json', 'code/q04/release_q04.py', 'code/q04/q04_expected.json', 'code/p6q05/release_q05.py', 'code/p6q05/q05_expected.json', 'code/p6q06/release_q06.py', 'code/p6q06/q06_expected.json', 'code/p6q06/q07/release_q07.py', 'code/p6q06/q07/q07_expected.json')}
        rep = {"landmark": "XVI-Q", "release": "Q07", "zip_sha256": hashlib.sha256(zp.read_bytes()).hexdigest(),
               "source_tree_sha256": SS.tree_sha(root), "tree_matches_reference": SS.tree_sha(root, exclude=(REF,))
               == ref["tree_sha256_without_reference"], "full_build_ok": built, "compared": len(ref["values"]),
               "equal": len(ref["values"]) - len(diff), "differences": diff, "gate_and_expected_sha256": v_files}
        rep["pass"] = bool(built and not diff and rep["tree_matches_reference"])
    out = Path(sys.argv[sys.argv.index("--report") + 1]) if "--report" in sys.argv else Q / "build" / "rebuild_report.json"
    out.write_text(json.dumps(rep, indent=1) + "\n")
    print(json.dumps(rep, indent=1))
    return 0 if rep["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
