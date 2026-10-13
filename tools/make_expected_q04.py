#!/usr/bin/env python3
"""AUTHOR tool (freeze): code/q04/q04_expected.json = pinned copies of every Q01..Q04 certificate ('seconds' removed)
+ sha256 of code files; pins its sha into release_q04.EXPECTED_SHA256."""
import hashlib, json, re, sys  # noqa: E401
from pathlib import Path

Q = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Q / "code" / "q04"))
import release_q04 as R  # noqa: E402
NAMES = ["Q01_intake", "Q01_O2", "Q01_guard", "Q01_release", "Q02_bulk_graphs", "Q02_D12_reconciliation",
         "Q02_bulk_identity", "Q02_O1", "Q02_D04_scrutiny", "Q02_release", "Q03_O1_record", "Q03_guard_ring",
         "Q03_release", "Q04_P6_inventory", "Q04_P6_identity", "Q04_P6_eval", "Q04_P6_status", "Q04_guard_scan"]
FILES = sorted({p.relative_to(Q).as_posix() for p in (Q / "code").rglob("*.py") if "__pycache__" not in p.parts
                and p.name != "release_q04.py"} | {"code/q02/q02_expected.json", "code/q03/q03_expected.json",
                                                    "intake/imports.json", "intake/errata_p23.json",
                                                    "reproducibility/independent_rebuild_report_q03.json"}
               | {p.relative_to(Q).as_posix() for p in (Q / "tools").glob("*.py") if p.name != "make_expected_q04.py"})
exp = {"certificates": {n: R.strip(json.loads((Q / "certificates" / f"{n}.json").read_text())) for n in NAMES},
       "file_sha256": {f: hashlib.sha256((Q / f).read_bytes()).hexdigest() for f in FILES}}
(Q / "code/q04/q04_expected.json").write_text(json.dumps(exp, indent=1, sort_keys=True) + "\n")
s = hashlib.sha256((Q / "code/q04/q04_expected.json").read_bytes()).hexdigest()
g = Q / "code/q04/release_q04.py"
g.write_text(re.sub(r'EXPECTED_SHA256 = "[^"]*"', f'EXPECTED_SHA256 = "{s}"', g.read_text()))
print(s)
