"""XVI-Q Q10 source set: Q09 set + Q10 reference values, expected, new-manuscript intake (tools/source_set.py is pinned by
the Q04 layer gate and left unchanged)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import source_set_q06_q07_q08_q09 as S0  # noqa: E402

GLOBS = tuple(S0.GLOBS) + ("reproducibility/reference_values_q10.json", "code/p6q06/q07/q08/q09/q10/q10_expected.json",
                               "intake/new_manuscripts_q10.json")
wrap = S0.wrap


def source_files(root):
    root = Path(root)
    out = set()
    for g in GLOBS:
        out |= {p.relative_to(root).as_posix() for p in root.glob(g) if p.is_file() and "__pycache__" not in p.parts}
    return sorted(out)


def tree_sha(root, exclude=()):
    import hashlib
    root = Path(root)
    lines = [f"{r}\t{hashlib.sha256((root / r).read_bytes()).hexdigest()}" for r in source_files(root) if r not in exclude]
    return hashlib.sha256(("\n".join(lines) + "\n").encode()).hexdigest()


def make_zip(out, root):
    import zipfile
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for r in source_files(root):
            zi = zipfile.ZipInfo(f"XVI-Q/{r}", date_time=(2026, 10, 11, 0, 0, 0))
            zi.compress_type = zipfile.ZIP_DEFLATED
            z.writestr(zi, (Path(root) / r).read_bytes())
