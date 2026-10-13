#!/usr/bin/env python3
"""AUTHOR tool (freeze): code/p6q06/q07/q08/q09/q10/q10_expected.json = pinned copies of every Q01..Q10 certificate ('seconds'
removed) + sha256 of EVERY code/**/*.py, tools/*.py, tests/*.py (except the Q10 gate itself) and the intake /
expected files; pins its sha into release_q10.EXPECTED_SHA256."""
import hashlib, json, re, sys  # noqa: E401
from pathlib import Path

Q = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Q / "code" / "p6q06" / "q07" / "q08" / "q09"))
sys.path.insert(0, str(Q / "code" / "p6q06" / "q07" / "q08" / "q09" / "q10"))
sys.path.insert(0, str(Q / "code" / "p6q06" / "q07" / "q08"))
sys.path.insert(0, str(Q / "code" / "p6q06" / "q07"))
sys.path.insert(0, str(Q / "code" / "p6q06"))
import release_q10 as R  # noqa: E402
import release_q07 as R7  # noqa: E402
import release_q06 as R6  # noqa: E402
NAMES = ["Q01_intake", "Q01_O2", "Q01_guard", "Q01_release", "Q02_bulk_graphs", "Q02_D12_reconciliation",
         "Q02_bulk_identity", "Q02_O1", "Q02_D04_scrutiny", "Q02_release", "Q03_O1_record", "Q03_guard_ring",
         "Q03_release", "Q04_P6_inventory", "Q04_P6_identity", "Q04_P6_eval", "Q04_P6_status", "Q04_guard_scan",
         "Q04_release", "Q05_rules_validation", "Q05_P6_graphs", "Q05_P6_rows", "Q05_P6_status", "Q05_guard_scan",
         "Q05_release"] + R6.Q06_CERTS + ["Q06_release"] + R.Q07_CERTS + R.Q08_CERTS + R.Q09_CERTS + ["Q09_release"] + R.Q10_CERTS
FILES = sorted(R.code_files() | {"code/q02/q02_expected.json", "code/q03/q03_expected.json", "code/q04/q04_expected.json",
                                 "code/p6q05/q05_expected.json", "code/p6q06/q06_expected.json", "code/p6q06/q07/q07_expected.json", "code/p6q06/q07/q08/q08_expected.json", "code/p6q06/q07/q08/q09/guard_scan3_rules.json", "code/p6q06/q07/q08/q09/q09_expected.json", "intake/new_manuscripts_q10.json", "intake/imports.json",
                                 "intake/errata_p23.json"})
exp = {"certificates": {n: R.strip(json.loads((Q / "certificates" / f"{n}.json").read_text())) for n in NAMES},
       "file_sha256": {f: hashlib.sha256((Q / f).read_bytes()).hexdigest() for f in FILES}}
p = Q / "code/p6q06/q07/q08/q09/q10/q10_expected.json"
p.write_text(json.dumps(exp, indent=1, sort_keys=True) + "\n")
s = hashlib.sha256(p.read_bytes()).hexdigest()
g = Q / "code/p6q06/q07/q08/q09/q10/release_q10.py"
g.write_text(re.sub(r'EXPECTED_SHA256 = "[^"]*"', f'EXPECTED_SHA256 = "{s}"', g.read_text()))
print(s)
