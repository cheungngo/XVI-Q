"""XVI-Q source set (what the source-only zip contains) and helpers."""
from __future__ import annotations
import hashlib, textwrap, zipfile  # noqa: E401
from pathlib import Path

GLOBS = ("README.md", "REBUILD_INSTRUCTIONS.md", "registration/*.json", "intake/*.json", "provenance/**/*",
         "code/**/*.py", "code/*/q0*_expected.json", "reproducibility/independent_rebuild_report_q02.json", "reproducibility/independent_rebuild_report_q03.json", "tests/*.py", "tools/*.py", "manuscript/notes/*.tex", "reproducibility/reference_values_q04.json")


def source_files(root):
    root = Path(root)
    out = set()
    for g in GLOBS:
        out |= {p.relative_to(root).as_posix() for p in root.glob(g) if p.is_file() and "__pycache__" not in p.parts}
    return sorted(out)


def tree_sha(root, exclude=()):
    root = Path(root)
    lines = [f"{r}\t{hashlib.sha256((root / r).read_bytes()).hexdigest()}" for r in source_files(root) if r not in exclude]
    return hashlib.sha256(("\n".join(lines) + "\n").encode()).hexdigest()


def make_zip(out, root):
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for r in source_files(root):
            zi = zipfile.ZipInfo(f"XVI-Q/{r}", date_time=(2026, 10, 10, 0, 0, 0))
            zi.compress_type = zipfile.ZIP_DEFLATED
            z.writestr(zi, (Path(root) / r).read_bytes())


def wrap(s, w=150):
    s = s.replace("\\end{verbatim}", "\\end {verbatim}")   # a listed file must not close the listing
    s = "".join(ch if ord(ch) < 128 else f"<U+{ord(ch):04X}>" for ch in s)   # Q02: pdflatex verbatim is ASCII-only
    return "\n".join(textwrap.fill(l, w, break_long_words=True, replace_whitespace=False) if len(l) > w else l
                     for l in s.splitlines())
