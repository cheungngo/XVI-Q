"""XVI-Q d guard (Q01; Q02 revision: generator-only Poly refusal, apart removed from the allow-list): XVI-P p23 dred (imported read-only, imports/XVI-P-p23) + ONE entry-point rule.

p23 review bypasses (sympy.polys.rem namespace, polytools.rem captured before import, reduced, groebner,
Matrix.subs d->Dummy then Dummy->4, Poly(...).rep.eval) all reach either the polynomial constructor Poly.new (class
attribute, looked up at call time -> namespace aliases and pre-import captures cannot avoid it) or an internal rename
of d.  Rule: such an operation touching d is admitted only if the ENTRY function -- the first sympy frame entered from
non-sympy code -- is an allow-listed rational-normalisation entry (cancel/factor/together/apart/simplify/expand/...
and arithmetic dunders), or inside dred.authorised_d.  Everything else raises PrematureD4Substitution (deny by default).
"""
from __future__ import annotations
import sys
from pathlib import Path

QROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(QROOT / "imports" / "XVI-P-p23" / "code" / "p5_entry"))
import sympy as sp  # noqa: E402
import sympy.polys.polytools as _pt  # noqa: E402
import dred  # noqa: E402   (XVI-P p23, pristine import)

d = dred.d
PrematureD4Substitution = dred.PrematureD4Substitution
ALLOW_ENTRY = frozenset({
    "cancel", "factor", "factor_list", "factor_terms", "together", "simplify", "expand",
    "expand_mul", "ratsimp", "radsimp", "gcd", "lcm", "gcd_terms", "numer", "denom", "fraction", "as_numer_denom",
    "nsimplify_none", "count_ops", "equals", "is_zero", "could_extract_minus_sign",
    "__add__", "__radd__", "__sub__", "__rsub__", "__mul__", "__rmul__", "__truediv__", "__rtruediv__", "__pow__",
    "__rpow__", "__neg__", "__eq__", "__ne__", "__hash__", "__new__", "_eval_power", "flatten",
})
_INTERNAL_MARK = (f"{Path(sp.__file__).parent}", "mpmath", str(Path(dred.__file__).resolve()), str(Path(__file__).resolve()))
ENTRY_POINTS = ["Poly.new (polynomial construction with d as generator or in the domain)",
                "dred._check_rule (any internal rename/substitution of d)"]


def _is_internal(fr):
    f = fr.f_code.co_filename
    return any(m in f for m in _INTERNAL_MARK)


WRAPPERS = frozenset({"threaded_func", "wrapper", "_sympifyit_wrapper", "__sympifyit_wrapper", "binary_op_wrapper",
                      "_func", "wrapped", "inner", "decorator", "cacheit", "wrapper_cache"})
POLY_CLASSES = frozenset({"Poly", "PurePoly", "DMP", "DMP_Python", "DUP_Flint", "PolyElement", "PolyRing"})


def _entry_frame():
    """first sympy frame entered from non-sympy code, skipping decorator-wrapper frames."""
    fr = sys._getframe(2)
    chain = []
    while fr is not None and _is_internal(fr):
        chain.append(fr)
        fr = fr.f_back
    for f in reversed(chain):                       # outermost internal frame first
        if f.f_code.co_name not in WRAPPERS:
            return f
    return chain[-1] if chain else None


def entry_function():
    f = _entry_frame()
    return getattr(f.f_code, "co_qualname", f.f_code.co_name) if f is not None else None


def _entry_allowed():
    q = entry_function()
    if q is None:
        return False
    if q.split(".")[0] in POLY_CLASSES:             # Poly(...), Poly.method(...) are never an allowed entry
        return False
    return q.split(".")[-1] in ALLOW_ENTRY


def _authorised():
    return bool(getattr(dred._AUTH, "active", False))


def _deny(what):
    raise PrematureD4Substitution(f"XVI-Q guard: {what} (entry {entry_function()!r} not allow-listed); use "
                                  "dred.substitute_d / dred.authorised_d")


def _touches_d(gens, dom):
    """Q02 revision: refuse only when d is a polynomial GENERATOR (every reviewer bypass -- rem/reduced/groebner/
    compose/resultant/shift/rep.eval/apart -- evaluates or divides in the d variable).  d inside the coefficient
    domain (e.g. Poly(expr, x) over QQ(d), as in XVI-P integrals.fp_two_point) cannot set d = 4 and is admitted."""
    return any(dred._has_d(g) for g in gens)


_orig_new = _pt.Poly.new.__func__


def _guarded_new(cls, rep, *gens):
    if not _authorised():
        dom = getattr(rep, "dom", None)
        if _touches_d(gens, dom) and not _entry_allowed():
            _deny("Poly construction with d in the polynomial domain")
    return _orig_new(cls, rep, *gens)


_pt.Poly.new = classmethod(_guarded_new)
_p23_rule = dred._check_rule


def _check_rule(pairs, what, depth=2):
    if _authorised():
        return
    for old, new in pairs:
        if dred._has_d(old) and not _entry_allowed():
            _deny(f"{what} of {old} -> {new}")
    _p23_rule(pairs, what, depth + 1)


dred._check_rule = _check_rule

# Q03: PolyRing / FracField with d as a GENERATOR (sympy.polys.rings.ring / fields.field / sring, then
# PolyElement.evaluate(d, 4)) bypassed Poly.new (Q02 independent review).  Constructing such rings stays admitted (QQ[d]
# is the legitimate coefficient domain of Poly(expr, x)); every way of assigning a VALUE to a generator of a ring/field
# that has d among its symbols (evaluate / subs / __call__ / compose, on PolyElement and FracElement) is refused unless
# the entry is allow-listed or inside dred.authorised_d.
import sympy.polys.rings as _rings  # noqa: E402
import sympy.polys.fields as _fields  # noqa: E402


def _ring_has_d(el):
    r = getattr(el, "ring", None) or getattr(el, "field", None)
    return r is not None and any(dred._has_d(sp.sympify(x)) for x in getattr(r, "symbols", ()))


EVAL_ENTRY = frozenset({"evaluate", "subs", "__call__", "compose"})   # internal sympy use (e.g. heugcd evaluation
# points inside an allow-listed cancel/apart) is admitted: it never returns a d = 4 value to the caller.


def _wrap_eval(klass, meth):
    orig = getattr(klass, meth)

    def _m(self, *a, **k):
        direct = not _is_internal(sys._getframe(1))        # called from non-sympy code (the reviewer bypass)
        if not _authorised() and _ring_has_d(self) and (direct or (entry_function() or "").split(".")[-1] in EVAL_ENTRY):
            _deny(f"{klass.__name__}.{meth} on a ring/field with d as a generator")
        return orig(self, *a, **k)
    _m.__name__ = meth
    setattr(klass, meth, _m)


for _k in (_rings.PolyElement, _fields.FracElement):
    for _meth in ("evaluate", "subs", "__call__", "compose"):
        if hasattr(_k, _meth):
            _wrap_eval(_k, _meth)
ENTRY_POINTS.append("PolyElement/FracElement evaluate/subs/__call__/compose on a ring/field with d as generator")
INSTALLED = True
