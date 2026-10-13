#!/usr/bin/env python3
"""XVI-Q Q05 hardened static policy scan (Q04 independent review: 7 evasion forms; scan limited to code/q0*/*.py).
Scans EVERY .py under code/ (recursive).  Flags, outside the declared probe files:
  F1 any reference (call, alias import, attribute, bare name, or string constant) to an escape primitive
     (srepr, parse_expr, lambdify, dup_eval*, dmp_eval*, as_coefficients_dict, evalf, N) -- closes aliased imports
     and getattr(x, 'srepr')
  F2 getattr / __getattribute__ with a NON-literal attribute name (string concatenation)
  F3 .replace(...) whose receiver is not a string literal, a Path, or allow-listed -- closes str/f-string/format/__str__
     rename-and-reparse
  F4 rebinding of builtins str / repr / format / getattr
  F5 name-based symbol identification: comparisons on `.name` (tree walk recognising d by name)
CLAIM (honest, downgraded): best-effort POLICY scan of the package's own code.  It is NOT a security boundary: Python
is dynamic, so determined code can still evade it (e.g. exec of a computed string, importlib with computed names);
F6 flags exec/eval/compile/importlib/__import__ to raise the cost, but the guarantee is review + scan, not proof."""
from __future__ import annotations
import ast, json, sys  # noqa: E401
from pathlib import Path

Q = Path(__file__).resolve().parents[2]
PRIM = {"srepr", "parse_expr", "lambdify", "dup_eval", "dmp_eval", "dmp_eval_in", "dmp_eval_tail", "dup_eval_in",
        "as_coefficients_dict", "evalf", "N"}
BUILTINS = {"str", "repr", "format", "getattr"}
DYN = {"exec", "eval", "compile", "__import__", "import_module"}
PROBE_FILES = {"code/q04/guard_scan.py", "code/q01/guard_probe.py", "code/p6q05/guard_scan2.py"}
# .replace receivers that are provably plain text manipulation of non-expression data (file text / labels)
REPLACE_OK = {}
# reviewed exceptions (line-pinned; the gate pins file sha256, so a moved/edited line re-triggers review):
REVIEWED = {("code/q01/o2_q01.py", 34): "__import__('zipfile') for provenance zip hashing; no expression data",
            ("code/q01/qguard.py", 132): "getattr(klass, meth) inside the guard's own wrapper installation"}
# 'N' as a bare name is a common local (numerators); sympy.N is caught via import alias / attribute (.N) instead.


def files():
    return sorted(p for p in (Q / "code").rglob("*.py") if "__pycache__" not in p.parts)


def _hits(rel, tree):
    h = []
    for n in ast.walk(tree):
        if isinstance(n, (ast.Name,)) and n.id in (PRIM - {"N"}) | DYN and (rel, n.lineno) not in REVIEWED:
            h.append(f"{rel}:{n.lineno}: F1 name {n.id}")
        if isinstance(n, ast.Attribute) and n.attr in PRIM | {"import_module"}:
            h.append(f"{rel}:{n.lineno}: F1 attr .{n.attr}")
        if isinstance(n, ast.alias) and (n.name.split(".")[-1] in PRIM | DYN):
            h.append(f"{rel}:{getattr(n, 'lineno', 0)}: F1 import {n.name}")
        if isinstance(n, ast.Constant) and isinstance(n.value, str) and n.value in PRIM:
            h.append(f"{rel}:{n.lineno}: F1 string '{n.value}'")
        if isinstance(n, ast.Call):
            f = n.func
            nm = f.attr if isinstance(f, ast.Attribute) else f.id if isinstance(f, ast.Name) else None
            if nm in ("getattr", "__getattribute__") and len(n.args) >= 2 and not isinstance(n.args[1], ast.Constant) \
                    and (rel, n.lineno) not in REVIEWED:
                h.append(f"{rel}:{n.lineno}: F2 getattr non-literal")
            if nm == "replace" and isinstance(f, ast.Attribute):
                recv = f.value
                lit = isinstance(recv, ast.Constant)
                if not lit and (rel, n.lineno) not in REPLACE_OK:
                    h.append(f"{rel}:{n.lineno}: F3 .replace on non-literal")
        if isinstance(n, (ast.Assign, ast.AnnAssign, ast.AugAssign)):
            tg = n.targets if isinstance(n, ast.Assign) else [n.target]
            for t in tg:
                if isinstance(t, ast.Name) and t.id in BUILTINS:
                    h.append(f"{rel}:{n.lineno}: F4 rebinding {t.id}")
        if isinstance(n, (ast.FunctionDef, ast.arg)) and getattr(n, "name", getattr(n, "arg", None)) in BUILTINS:
            h.append(f"{rel}:{n.lineno}: F4 rebinding {getattr(n, 'name', getattr(n, 'arg', None))}")
        if isinstance(n, ast.Compare):
            for side in [n.left] + n.comparators:
                if isinstance(side, ast.Attribute) and side.attr == "name":
                    h.append(f"{rel}:{n.lineno}: F5 compare on .name")
    return h


def scan():
    out = []
    for p in files():
        rel = p.relative_to(Q).as_posix()
        if rel in PROBE_FILES:
            continue
        out += _hits(rel, ast.parse(p.read_text()))
    return out


EVASIONS = {   # the reviewer's 7 forms, each must be flagged
    "aliased_lambdify": "from sympy import lambdify as f\ng = f(x, e)\n",
    "getattr_concat_srepr": "g = getattr(sp, 'sr' + 'epr')(e)\n",
    "fstring_replace": "v = sp.sympify(f'{e}'.replace('d', '4'))\n",
    "dunder_str_replace": "v = sp.sympify(e.__str__().replace('d', '4'))\n",
    "format_replace": "v = sp.sympify(format(e).replace('d', '4'))\n",
    "rebound_str": "str = lambda x: x\nv = str(e)\n",
    "name_based_walk": "def w(x):\n    if x.name == 'd':\n        return 4\n",
}


def selftest():
    return {k: bool(_hits("<probe>", ast.parse(src))) for k, src in EVASIONS.items()}


def report():
    h = scan()
    st = selftest()
    return {"scope": "every .py under code/ (recursive), minus declared probe files", "probe_files": sorted(PROBE_FILES),
            "n_files_scanned": len(files()) - len([p for p in files() if p.relative_to(Q).as_posix() in PROBE_FILES]),
            "rules": ["F1 primitive reference/alias/string", "F2 non-literal getattr", "F3 .replace on non-literal",
                      "F4 builtin rebinding", "F5 .name comparison", "F6 dynamic exec/eval/import"],
            "reviewer_evasions_flagged": st, "reviewed_exceptions": {f"{a}:{b}": v for (a, b), v in REVIEWED.items()}, "all_reviewer_evasions_flagged": all(st.values()),
            "scan_hits": h, "scan_clean": not h,
            "claim": "best-effort policy scan; NOT a security boundary (dynamic Python can still evade it)",
            "runtime_closable": False}


if __name__ == "__main__":
    r = report()
    print(json.dumps({k: r[k] for k in ("n_files_scanned", "all_reviewer_evasions_flagged", "scan_clean")}), *r["scan_hits"], sep="\n")
