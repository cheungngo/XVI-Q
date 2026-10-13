#!/usr/bin/env python3
"""XVI-Q Q07 goal 1 (main test): do LOCAL evanescent kappa^3 counterterms exist that make the ghost rows vanish
identically?  Inputs: kappa^2 two-point poles with the p22 evanescent CTs applied (4D-projected Gamma_{Kc},
Gamma_{cbar c}); kappa^3 poles Gamma_{X c h} (X = K / cbar) and Gamma_{c* c c} plus a GENERIC local evanescent
counterterm each: every monomial of the right polarisation multidegree and momentum degree containing at least one T
(tilde dot) between polarisations, with a free coefficient.  The row is linear, so existence = solvability of a linear
system over Q(i).  Solvable -> the evanescent remainder is sigma0-exact (cancellable); not solvable -> explicit
evanescent obstruction on these rows (residual of the least-squares-free reduction recorded)."""
from __future__ import annotations
import itertools, sys  # noqa: E401
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "p6q05"))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import rows as R  # noqa: E402
import rows6  # noqa: E402
import uvpole as U  # noqa: E402
from evct import four_d  # noqa: E402

MOM = ("p1", "p2")


def monomials(pols, nmom):
    """all products of H/T dot symbols with exact degree pols[name] in each polarisation and total momentum degree
    nmom, at least one T (T only between polarisations)."""
    names = sorted(pols)
    slots = []
    for n in names:
        slots += [n] * pols[n]
    out = set()

    def rec(rem, mom_left, acc, hasT):
        if not rem:
            if mom_left % 2 == 0 and hasT:
                k = mom_left // 2
                for pp in itertools.combinations_with_replacement([U.H("p1", "p1"), U.H("p1", "p2"), U.H("p2", "p2")], k):
                    out.add(acc * sp.Mul(*pp))
            return
        x = rem[0]
        rest = list(rem[1:])
        # pair x with a momentum
        if mom_left > 0:
            for m in MOM:
                rec(rest, mom_left - 1, acc * U.H(x, m), hasT)
        # pair x with another remaining slot
        seen = set()
        for j, y in enumerate(rest):
            if y in seen:
                continue
            seen.add(y)
            r2 = rest[:j] + rest[j + 1:]
            rec(r2, mom_left, acc * U.H(x, y), hasT)
            rec(r2, mom_left, acc * U.T(x, y), True)
    rec(slots, nmom, sp.Integer(1), False)
    return sorted(out, key=str)


def ansatz(pols, nmom, tag):
    ms = monomials(pols, nmom)
    xs = sp.symbols(f"{tag}0:{len(ms)}")
    return sp.expand(sp.I * sum(x * m for x, m in zip(xs, ms))), list(xs)


def run(target):
    if target == "ccK":
        A, xa = ansatz({"aK": 2, "gc1": 1, "a3": 2}, 3, "xa")
        loopX = ("G03", "G07")
    else:
        A, xa = ansatz({"gb2": 1, "gc1": 1, "a3": 2}, 4, "xa")
        loopX = ("G04", "G09")
    B, xb = ansatz({"z": 1, "gc1": 1, "gc2": 1}, 3, "xb")
    orig = (R.load, R.D13_own, R.D06_own, R.D13_K_ext)
    try:
        def load(g):
            e = orig[0](g)
            if g == loopX[0]:
                return sp.expand(e + A)
            if g == "G10":
                return sp.expand(e + B)
            return e
        R.load = load
        R.D13_own = lambda *a: four_d(orig[1](*a))                  # noqa: E731
        R.D06_own = lambda *a: four_d(orig[2](*a))                  # noqa: E731
        R.D13_K_ext = lambda *a: four_d(orig[3](*a))                # noqa: E731
        P = rows6.parts(target)
    finally:
        R.load, R.D13_own, R.D06_own, R.D13_K_ext = orig
    C = rows6.coeffs(target)
    tot = sp.expand(sum(C[k] * P[k] for k in P))
    xs = xa + xb
    gens = sorted([s for s in tot.free_symbols if s not in set(xs)], key=str)
    eqs = sp.Poly(tot, *gens).coeffs()
    sol = sp.linsolve(eqs, xs)
    ok = sol != sp.EmptySet
    res = {"row": target, "n_ct_params": len(xs), "n_equations": len(eqs), "solvable": ok}
    if ok:
        s0 = list(sol)[0]
        free = set().union(*[sp.sympify(v).free_symbols for v in s0]) & set(xs)
        res["n_free_params"] = len(free)
        part = {x: (sp.Integer(0) if x in free else sp.sympify(v)) for x, v in zip(xs, s0)}
        ctA = sp.expand(A.xreplace({x: v.xreplace({f: 0 for f in free}) for x, v in part.items() if x in xa}))
        ctB = sp.expand(B.xreplace({x: v.xreplace({f: 0 for f in free}) for x, v in part.items() if x in xb}))
        res["particular_CT_Gamma_Xch"] = str(ctA)
        res["particular_CT_Gamma_cstar_cc"] = str(ctB)
        res["particular_CT_Xch_n_terms"] = len(sp.Add.make_args(ctA)) if ctA != 0 else 0
        res["particular_CT_cstar_cc_n_terms"] = len(sp.Add.make_args(ctB)) if ctB != 0 else 0
    else:
        # rank deficiency report: solve with the remainder only on the evanescent part
        res["note"] = "no local evanescent kappa^3 CT in the ansatz removes the remainder"
    return res


if __name__ == "__main__":
    for t in ("ccK", "ccb"):
        r = run(t)
        print({k: (v if not isinstance(v, str) or len(v) < 200 else v[:200]) for k, v in r.items()}, flush=True)
