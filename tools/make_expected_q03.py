#!/usr/bin/env python3
"""AUTHOR tool (run at freeze, after make_expected_q02): write code/q03/q03_expected.json = pinned copies of every
Q01/Q02/Q03 certificate ('seconds' removed) + sha256 of code/conventions files; pin its sha into release_q03."""
import hashlib, json, re, sys  # noqa: E401
from pathlib import Path

Q = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Q / "code" / "q03"))
import release_q03 as R  # noqa: E402
NAMES = ["Q01_intake", "Q01_O2", "Q01_guard", "Q01_release", "Q02_bulk_graphs", "Q02_D12_reconciliation",
         "Q02_bulk_identity", "Q02_O1", "Q02_D04_scrutiny", "Q02_release", "Q03_O1_record", "Q03_guard_ring"]
FILES = ["code/q01/q01_conventions.py", "code/q01/release_q01.py", "code/q01/qguard.py", "code/q02/bulk_graphs.py",
         "code/q02/bulk_identity.py", "code/q02/d04_route4.py", "code/q02/release_q02.py", "code/q02/q02_expected.json",
         "code/q03/q03_certs.py", "intake/imports.json", "intake/errata_p23.json", "registration/local_registration.json"]
exp = {"certificates": {n: R.strip(json.loads((Q / "certificates" / f"{n}.json").read_text())) for n in NAMES},
       "file_sha256": {f: hashlib.sha256((Q / f).read_bytes()).hexdigest() for f in FILES}}
(Q / "code/q03/q03_expected.json").write_text(json.dumps(exp, indent=1, sort_keys=True) + "\n")
s = hashlib.sha256((Q / "code/q03/q03_expected.json").read_bytes()).hexdigest()
g = Q / "code/q03/release_q03.py"
g.write_text(re.sub(r'EXPECTED_SHA256 = "[^"]*"', f'EXPECTED_SHA256 = "{s}"', g.read_text()))
print(s)
