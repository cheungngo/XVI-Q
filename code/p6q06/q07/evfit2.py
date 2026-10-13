#!/usr/bin/env python3
"""XVI-Q Q07 goals 1+3: joint test whether LOCAL evanescent kappa^3 counterterms X exist with sigma0 X = - (evanescent
remainder) on ALL kappa^3 three-point rows at once (c.c.K, c.c.cbar, h.h.c; c.phi.phi has no remainder).
CT = -Gamma_pole + X.  Unknowns: generic evanescent (>= one T) local monomials for Gamma_{Kch}, Gamma_{cbar ch},
Gamma_{c* cc}, Gamma_hhh with the right multidegree; the SAME Gamma_{Kch} / Gamma_{cbar ch} counterterm enters the
ghost rows and the h.h.c row.  kappa^2 poles enter unmodified (their CT = -pole, p22, is sigma0-closed at kappa^2).
Linear system over Q(i): solvable -> the kappa^3 three-point defect is cancellable (pole level); otherwise the
reduced residual is an explicit evanescent obstruction on these rows."""
from __future__ import annotations
import json, sys, time  # noqa: E401
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
Q = HERE.parents[2]
sys.path.insert(0, str(HERE.parents[1] / "p6q05"))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import rows as R  # noqa: E402
import rows6  # noqa: E402
import gravity as G  # noqa: E402
import hhc  # noqa: E402
from evfit import ansatz  # noqa: E402
LOG = Q / "build" / "logs"


def ghost_row_total(target, XA, XB):
    loopX = ("G03", "G07") if target == "ccK" else ("G04", "G09")
    orig = R.load
    try:
        def load(g):
            e = orig(g)
            if g == loopX[0]:
                return sp.expand(e + XA)
            if g == "G10":
                return sp.expand(e + XB)
            return e
        R.load = load
        P = rows6.parts(target)
    finally:
        R.load = orig
    C = rows6.coeffs(target)
    return sp.expand(sum(C[k] * P[k] for k in P))


def hhc_ansatz_total(XH, XK, XC, sig):
    """row h.h.c contributions of the ansatz CTs (channels A, D, E) + the stored loop residual."""
    base = sp.sympify(json.loads((LOG / "q07_hhc_row.json").read_text())["residual"])
    A = R.strip_t(-sp.I * R.polarise(sp.expand(XH), "a1", "p1", "gc1"))
    D = 0
    E = 0
    for i in ("2", "3"):
        ai, ki, Fi = hhc.EXT[i]
        Gk = R.polarise(sp.expand(hhc._relabel(XK, i)), "aK", "KU", "KW")
        eh, gf = G.S2_EH_gf(hhc.ext(ai), G.probe_mat("E"), ki)
        D += R.strip_t(R.contract_id(Gk, sp.I * (eh + gf), "K", "E", "h"))
        E += R.strip_t(hhc.subst_vec(sp.expand(hhc._relabel(XC, i)), "gb2", G.Dvec(hhc.ext(ai), Fi)))
    return sp.expand(base + A + D + sig * E)


def main():
    t0 = time.time()
    XK, xk = ansatz({"aK": 2, "gc1": 1, "a3": 2}, 3, "xk")
    XC, xc = ansatz({"gb2": 1, "gc1": 1, "a3": 2}, 4, "xc")
    XS, xs = ansatz({"z": 1, "gc1": 1, "gc2": 1}, 3, "xs")
    XH, xh = ansatz({"a1": 2, "a2": 2, "a3": 2}, 4, "xh")
    sig = sp.sympify(json.loads((LOG / "q07_hhc_tree.json").read_text())["sigma_F"])
    rows = {"ccK": ghost_row_total("ccK", XK, XS), "ccb": ghost_row_total("ccb", XC, XS),
            "hhc": hhc_ansatz_total(XH, XK, XC, sig)}
    xall = xk + xc + xs + xh
    X = set(xall)
    eqs = []
    for r in rows.values():
        gens = sorted([s for s in r.free_symbols if s not in X], key=str)
        eqs += sp.Poly(r, *gens).coeffs()
    # Bose symmetry of the local hhh counterterm (it must come from a local functional): (12) and (13) swaps
    for mp in ({"p1": R.P2, "p2": R.P1, "a1": "a2", "a2": "a1"}, {"p1": R.P3, "p2": R.P2, "a1": "a3", "a3": "a1"}):
        d = sp.expand(XH - R.remap(XH, mp))
        if d != 0:
            gens = sorted([s_ for s_ in d.free_symbols if s_ not in X], key=str)
            eqs += sp.Poly(d, *gens).coeffs()
    sol = sp.linsolve(eqs, xall)
    ok = sol != sp.EmptySet
    out = {"unknowns": {"Gamma_Kch": len(xk), "Gamma_cbar_ch": len(xc), "Gamma_cstar_cc": len(xs), "Gamma_hhh": len(xh)},
           "n_equations": len(eqs), "solvable": ok, "seconds": round(time.time() - t0, 1)}
    if ok:
        s0 = list(sol)[0]
        free = set().union(*[sp.sympify(v).free_symbols for v in s0]) & X
        out["n_free"] = len(free)
        val = {x: sp.sympify(v).xreplace({f: 0 for f in free}) for x, v in zip(xall, s0)}
        for name, (Xe, xx) in {"Gamma_Kch": (XK, xk), "Gamma_cbar_ch": (XC, xc), "Gamma_cstar_cc": (XS, xs),
                               "Gamma_hhh": (XH, xh)}.items():
            e = sp.expand(Xe.xreplace({x: val[x] for x in xx}))
            out[f"particular_CT_{name}"] = str(e)
            out[f"particular_CT_{name}_n_terms"] = len(sp.Add.make_args(e)) if e != 0 else 0
        # verification: substitute the particular solution and recompute every row
        e = sp.expand(XH.xreplace({x: val[x] for x in xh}))
        out["hhh_CT_bose_symmetric"] = all(sp.expand(e - R.remap(e, mp)) == 0 for mp in (
            {"p1": R.P2, "p2": R.P1, "a1": "a2", "a2": "a1"}, {"p1": R.P3, "p2": R.P2, "a1": "a3", "a3": "a1"}))
        out["all_CT_evanescent"] = all(all(any(y.name.startswith("T_") for y in t.free_symbols)
                                           for t in sp.Add.make_args(sp.sympify(out[f"particular_CT_{n}"])))
                                       for n in ("Gamma_Kch", "Gamma_cbar_ch", "Gamma_hhh") if out[f"particular_CT_{n}"] != "0")
        out["verified_rows_zero"] = all(sp.expand(r.xreplace(val)) == 0 for r in rows.values())
    (LOG / "q07_evfit_joint.json").write_text(json.dumps(out, indent=1) + "\n")
    print({k: v for k, v in out.items() if not k.startswith("particular_CT_") or k.endswith("n_terms")})


if __name__ == "__main__":
    main()
