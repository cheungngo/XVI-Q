#!/usr/bin/env python3
"""XVI-Q Q04 d-guard: the 'leave sympy symbolics' bypass class (Q03 review finding guard_rename_class, low).

Runtime probes (recorded, NOT closable at runtime): once an expression is turned into Python data -- a pure-Python walk
over its terms with d replaced by 4, a str()/srepr() of it renamed and re-parsed, or raw dense coefficient lists
(densetools dup_eval/dmp_eval_in) -- d's value can be computed by code that never touches a sympy entry point.  No
runtime guard can see that.  It is therefore closed by POLICY + STATIC SCAN of the package's own code (code/q0*/*.py):
none of the escape primitives may appear there.  The gate requires the scan to be clean."""
from __future__ import annotations
import ast, sys  # noqa: E401
from pathlib import Path

Q = Path(__file__).resolve().parents[2]
FORBIDDEN_CALLS = {"srepr", "parse_expr", "dup_eval", "dmp_eval", "dmp_eval_in", "dmp_eval_tail", "dup_eval_in",
                   "as_coefficients_dict", "lambdify", "evalf", "N"}
SELF = {"code/q04/guard_scan.py", "code/q01/guard_probe.py"}   # probe files: call the primitives to show they raise / bypass


def scan():
    hits = []
    for p in sorted((Q / "code").glob("q0*/*.py")):
        rel = p.relative_to(Q).as_posix()
        if rel in SELF:
            continue
        tree = ast.parse(p.read_text())
        for n in ast.walk(tree):
            if isinstance(n, ast.Call):
                f = n.func
                name = f.attr if isinstance(f, ast.Attribute) else f.id if isinstance(f, ast.Name) else None
                if name in FORBIDDEN_CALLS:
                    hits.append(f"{rel}:{n.lineno}: {name}()")
                # str(...).replace(...) / repr(...).replace(...)  and  sympify(str(...)) / S(str(...))
                if name == "replace" and isinstance(f, ast.Attribute) and isinstance(f.value, ast.Call) and \
                        isinstance(f.value.func, ast.Name) and f.value.func.id in ("str", "repr"):
                    hits.append(f"{rel}:{n.lineno}: str(...).replace")
                if name in ("sympify", "S") and n.args and isinstance(n.args[0], ast.Call) and \
                        isinstance(n.args[0].func, ast.Name) and n.args[0].func.id in ("str", "repr"):
                    hits.append(f"{rel}:{n.lineno}: sympify(str(...))")
    return hits


def probes():
    """demonstrate the class on a probe expression (honest disclosure: these DO bypass the runtime guard)."""
    sys.path.insert(0, str(Q / "code" / "q01"))
    import qguard as G
    import sympy as sp
    d = G.d
    e = (d ** 2 - 16) / (d - 4) + 3 * d
    num, den = sp.fraction(sp.together(e))
    from fractions import Fraction

    def walk(x):                                     # pure-Python evaluation (Fractions), d -> 4: never re-enters sympy
        if isinstance(x, sp.Symbol) and x.name == d.name:     # identify d by NAME ('x == d' itself is refused by p22)
            return Fraction(4)
        if x.is_Integer or x.is_Rational:
            return Fraction(int(x.p), int(x.q))
        vals = [walk(a) for a in x.args]
        if x.is_Add:
            return sum(vals, Fraction(0))
        if x.is_Mul:
            out = Fraction(1)
            for v in vals:
                out *= v
            return out
        if x.is_Pow:
            return vals[0] ** int(vals[1])
        raise TypeError(type(x))
    res = {}
    try:
        _ = (d ** 2).args[0] == d
        res["python_term_walk_eq_d"] = "BYPASS"
    except G.PrematureD4Substitution:
        res["python_term_walk_eq_d"] = "raises (p22: equality test against d refused)"
    try:
        res["python_term_walk"] = f"BYPASS (value {walk(num)})"
    except G.PrematureD4Substitution:
        res["python_term_walk"] = "raises"
    try:
        v = sp.sympify(str(num).replace(str(d), "4"))
        res["str_rename_reparse"] = f"BYPASS (value {v})"
    except G.PrematureD4Substitution:
        res["str_rename_reparse"] = "raises"
    res["dense_list_eval"] = "BYPASS by construction (lists carry no symbol); disclosed since Q03"
    return res


def report():
    h = scan()
    return {"class": "leave-sympy (pure-Python term walk, str/srepr rename + re-parse, raw dense lists)",
            "runtime_closable": False,
            "why": "the computation happens outside every sympy entry point; a runtime guard cannot observe it",
            "probes": probes(), "policy": "package code (code/q0*/*.py) must not use any escape primitive",
            "forbidden": sorted(FORBIDDEN_CALLS) + ["str(...).replace", "sympify(str(...))"],
            "scan_hits": h, "scan_clean": not h,
            "n_files_scanned": len([p for p in (Q / "code").glob("q0*/*.py")]) - len(SELF)}


if __name__ == "__main__":
    import json
    print(json.dumps(report(), indent=1))
