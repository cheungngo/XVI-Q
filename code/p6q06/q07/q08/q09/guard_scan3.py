#!/usr/bin/env python3
"""XVI-Q Q09 scan extension (Q07 review, low): on top of p6q05/guard_scan2 (unchanged), flag
 (i) operator.methodcaller in any form (import, alias, attribute), and
 (ii) coefficient / exponent walks aimed at the opaque d: X.coeff(d..), X.as_coeff_exponent(d), X.as_independent(d),
     X.collect(d), X.as_poly(d), degree(X, d), and dict walks X.as_powers_dict()[d] / .get(d) /
     X.as_coefficients_dict()[d] / .get(d), where the argument is a d-like name (d, D, dim, dsym, d_sym, DIM, D_SYM)
     or an attribute ending in such a name.
Static and best-effort (the runtime guard stays primary); scope: every .py under code/, tools/, tests/ minus the
declared probe files of guard_scan2 and this file's probes."""
from __future__ import annotations
import ast, json, sys  # noqa: E401
from pathlib import Path

HERE = Path(__file__).resolve().parent
Q = HERE.parents[4]
sys.path.insert(0, str(HERE.parents[3] / "p6q05"))
import guard_scan2 as G2  # noqa: E402
DNAMES = {"d", "D", "dim", "dsym", "d_sym", "DIM", "D_SYM"}
WALK = {"coeff", "as_coeff_exponent", "as_independent", "collect", "as_poly", "degree"}
DICTS = set(json.loads((HERE / "guard_scan3_rules.json").read_text())["dict_walk_methods"])
SELF = "code/p6q06/q07/q08/q09/guard_scan3.py"


def _dlike(n):
    return (isinstance(n, ast.Name) and n.id in DNAMES) or (isinstance(n, ast.Attribute) and n.attr in DNAMES)


def _dictcall(n):
    return isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute) and n.func.attr in DICTS


def hits(tree):
    out = []
    for n in ast.walk(tree):
        if isinstance(n, ast.Name) and n.id == "methodcaller":
            out.append(("methodcaller", n.lineno))
        if isinstance(n, ast.Attribute) and n.attr == "methodcaller":
            out.append(("methodcaller", n.lineno))
        if isinstance(n, (ast.Import, ast.ImportFrom)) and {"methodcaller"} & {a.name for a in n.names}:
            out.append(("methodcaller_import", n.lineno))
        if isinstance(n, ast.Constant) and n.value == "methodcaller":
            out.append(("methodcaller_string", n.lineno))
        if isinstance(n, ast.Call):
            f = n.func
            name = f.attr if isinstance(f, ast.Attribute) else (f.id if isinstance(f, ast.Name) else None)
            if name in WALK and any(_dlike(a) for a in n.args):
                out.append((f"walk_{name}_d", n.lineno))
            if name == "get" and isinstance(f, ast.Attribute) and _dictcall(f.value) and any(_dlike(a) for a in n.args):
                out.append(("dict_walk_get_d", n.lineno))
        if isinstance(n, ast.Subscript) and _dictcall(n.value) and _dlike(n.slice):
            out.append(("dict_walk_subscript_d", n.lineno))
    return out


def files():
    fs = []
    for top in ("code", "tools", "tests"):
        fs += [p for p in (Q / top).rglob("*.py") if "__pycache__" not in p.parts]
    return sorted(p.relative_to(Q).as_posix() for p in fs if p.relative_to(Q).as_posix() not in (SELF,)
                  and p.relative_to(Q).as_posix() not in G2.PROBE_FILES)


PROBES = {
    "methodcaller_attr": "import operator\noperator.methodcaller('subs', d, 4)(e)",
    "methodcaller_alias": "from operator import methodcaller as mc\nmc('subs', d, 4)(e)",
    "methodcaller_getattr": "import operator\ngetattr(operator, 'methodcaller')('subs', d, 4)(e)",
    "coeff_d": "e.coeff(d, 1)",
    "coeff_exponent_d": "e.as_coeff_exponent(d)",
    "independent_d": "e.as_independent(d)",
    "degree_d": "sp.degree(e, d)",
    "powers_dict_subscript": "e.as_powers_dict()[d]",
    "powers_dict_get": "e.as_powers_dict().get(d)",
    "coefficients_dict_get": "e.as_coefficients_dict().get(qg.dsym)",
}


def report():
    hit = {}
    for rel in files():
        h = hits(ast.parse((Q / rel).read_text()))
        if h:
            hit[rel] = h
    probes = {k: bool(hits(ast.parse(v))) for k, v in PROBES.items()}
    base = G2.report()
    return {"scope": "every .py under code/, tools/, tests/ (recursive) minus declared probe files",
            "n_files_scanned": len(files()), "extension_hits": hit, "extension_clean": not hit,
            "probes_flagged": probes, "all_probes_flagged": all(probes.values()),
            "base_scan2_clean": base["scan_clean"], "base_scan2_evasions_flagged": base["all_reviewer_evasions_flagged"]}


if __name__ == "__main__":
    print(json.dumps(report(), indent=1))
