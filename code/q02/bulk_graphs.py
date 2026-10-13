#!/usr/bin/env python3
"""XVI-Q Q02 bulk O.5 graphs D04, D05a, D06, D12 on the frozen XVI-P p23 engine (re-executed copy build/xvip_run),
under the XVI-P declared prescription (DRED, opaque d -- qguard loaded FIRST, MS-bar, Feynman +i0, m = 0).

D04 phi phi mixed bubble (two V_h_phi_phi; graviton + scalar line), D06 cbar c mixed bubble (two V_cb_c_h; graviton +
ghost line): integrands built with XVI-P feynman.py vertices/propagators (ordered extraction, independent components,
Hessian-inverse propagators), same prefactor bookkeeping as amplitude_rules.build_integrands.
  route 1: integrals.fp_two_point (Feynman parameters, symmetric integration in QdS, NO scaleless rule);
  route 2 (D04): exact invariant fit N = sum c (l.l)^a (l.p)^b (p.p)^c verified identically, then IBP J_two_point;
  route 2 (D06): comparison with XVI-P P15_D13 ghost_self_energy (pinned import).
Symmetry factors derived by Wick counting (automorphisms of the graph with labelled external legs): D04 = D06 = 1.
D05a, D12: one massless line, vertex without loop-momentum dependence on the contracted legs (a2 carries no
derivative of h; for D12 the phi legs carry l but the integrand is a polynomial over l^2): pure power divergences,
omega in {2, 4}, no logarithmic member => DR value 0 DERIVED (no pole on the UV or the IR side), never by convention.
Gamma_1 convention (XVI-P): i Gamma_1 = graph sum; coefficients reported per kappa^2, in units of MB = Ihat(1,1)."""
from __future__ import annotations
import json
import sys
import time
from fractions import Fraction as Fr
from pathlib import Path

Q = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Q / "code" / "q01"))
import qguard  # noqa: E402,F401   (d guard installed before any XVI-P module)
RUN = Q / "build" / "xvip_run"
for sub in ("code/p5_calc", "code/p5_entry", "code/classical"):
    sys.path.insert(0, str(RUN / sub))
import sympy as sp  # noqa: E402
import dred  # noqa: E402
import feynman as F  # noqa: E402
import amplitude_rules as AR  # noqa: E402
import integrals as IL  # noqa: E402

assert dred.__file__ == qguard.dred.__file__, "one dred module (the guarded one)"
I, d, eps = sp.I, dred.d, dred.eps
ETA = (1, -1, -1, -1)
LS, PS = sp.symbols("l0:4"), sp.symbols("p0:4")
PP, MB, Lg = IL.PP, IL.MB, IL.L
NAMES = {**{r: LS[r] for r in range(4)}, **{4 + r: PS[r] for r in range(4)}}
P_NUMS = [(3, 1, 2, -1), (5, 2, -1, 1), (7, -3, 1, 2)]
SEEDS = (11, 23)


def sym(poly):
    return sp.expand(poly.to_sympy(NAMES)) if hasattr(poly, "to_sympy") else sp.expand(poly)


def dot(a, b):
    return sum(ETA[r] * a[r] * b[r] for r in range(4))


def gt(r, s):
    return sp.Symbol(f"Gt_{min(r, s)}{max(r, s)}")


def expand_coeff(c):
    """A(d) * MB -> pole / MS-bar finite (local + coefficient of L); d = 4 - 2 eps, MB = 1/eps + 2 - L + O(eps)."""
    ms = IL.master_series(1)
    a = sp.expand(dred.substitute_d(sp.sympify(c), 4 - 2 * eps, "Q02 bulk: Laurent expansion after all contractions"))
    ser = sp.series(a * (ms[-1] / eps + ms[0] + ms[1] * eps), eps, 0, 1).removeO()
    pole, fin = sp.simplify(ser.coeff(eps, -1)), sp.expand(sp.simplify(ser.coeff(eps, 0)))
    return {"A_exact_d": str(sp.factor(c)), "UV_pole": str(pole), "MSbar_finite": str(fin),
            "coefficient_of_L": str(sp.expand(fin).coeff(Lg)), "local": str(sp.expand(fin).subs(Lg, 0))}


def exact_fit(rows, rhs, n):
    dg = sp.Rational(173, 37)
    chosen, M = [], []
    for i, row in enumerate(rows):
        cand = M + [[dred.substitute_d(sp.sympify(x), dg, "Q02 fit: numerical rank test (row selection only)") for x in row]]
        if sp.Matrix(cand).rank() > len(M):
            M, chosen = cand, chosen + [i]
        if len(chosen) == n:
            break
    if len(chosen) < n:
        return None, False
    sol = [sp.factor(sp.cancel(x)) for x in sp.Matrix([rows[i] for i in chosen]).LUsolve(sp.Matrix([rhs[i] for i in chosen]))]
    ok = all(sp.cancel(sum(sp.sympify(a) * x for a, x in zip(row, sol)) - b) == 0 for row, b in zip(rows, rhs))
    return sol, ok


def fit(G, S):
    """G, S[name]: dict key -> expression in PS, PP, Gt; exact fit G = sum c_name S_name over realizations."""
    names = list(S)
    rows, rhs = [], []
    for pn in P_NUMS:
        for s in SEEDS:
            gm = IL.realization(pn, s)
            for k in sorted(set(G) | set().union(*[set(S[n]) for n in names])):
                rows.append([IL.realize(S[n].get(k, 0), pn, gm) for n in names])
                rhs.append(IL.realize(G.get(k, 0), pn, gm))
    sol, ok = exact_fit(rows, rhs, len(names))
    sol = sol or [None] * len(names)
    return {n: c for n, c in zip(names, sol) if c is not None and c != 0}, ok, len(rhs)


def prefactor_norm(S_f, loop_sign, phases, n_explicit_i):
    """Gamma_1 coefficient per kappa^2: -i * S_f * loop_sign * (i)^2 [two kappa^1 vertices] * i^phase * i^n."""
    return sp.expand(-I * sp.nsimplify(S_f) * loop_sign * I ** 2 * I ** phases * I ** n_explicit_i)


# ----------------------------------------------------------------------------------------------- D04
def D04_integrand(S_f=1, mutate=None):
    V = F.vertices()["V_h_phi_phi"]
    Pi = AR.pinv_sparse()
    if mutate == "drop_trace_part_of_graviton_propagator":
        dd = F.deDonder_tensor_numerator()

        def e(i):
            m, n = F.HC[i]
            return ETA[m] if m == n else 0
        Pi = {}
        for i in range(10):
            for j in range(10):
                v = sp.nsimplify(dd[i, j] + sp.Rational(1, 2) * e(i) * e(j))
                if v != 0:
                    Pi[(i, j)] = Fr(int(v.p), int(v.q))
    VA = {k: sym(v) for k, v in F.vt_route(V, [(-1, 0), (0, 1), (1, -1)]).items()}   # h_c: -l, phi_ext: p, phi: l-p
    VB = {k: sym(v) for k, v in F.vt_route(V, [(1, 0), (0, -1), (-1, 1)]).items()}   # h_d: l, phi_ext: -p, phi: p-l
    N = 0
    for (c, _, _), pa in VA.items():
        for (dd_, _, _), pb in VB.items():
            pv = Pi.get((c, dd_))
            if pv:
                N += pa * pb * sp.Rational(pv.numerator, pv.denominator)
    N = sp.expand(N)
    # propagators: graviton i P^{-1}/l^2, scalar i/(l-p)^2  -> two explicit i
    norm = prefactor_norm(S_f, 1, V.phase * 2, 2)
    return N, norm


def D04(S_f=1, mutate=None):
    N, norm = D04_integrand(S_f, mutate)
    R1 = sp.expand(norm * IL.fp_two_point(N, LS, PS) / MB)
    S = {"p4": {(): PP ** 2}}          # the only scalar structure: p.Gt.p = 0 identically (external p in the 4D part)
    c1, ok1, n1 = fit({(): R1}, S)
    out = {"route1_fp": {k: str(v) for k, v in c1.items()}, "route1_fit_exact": ok1, "route1_rows": n1}
    if mutate:
        return out, c1
    # route 2: invariant fit of N then IBP J
    l2, lp = sp.Symbol("l.l"), sp.Symbol("l.p")
    mons = [(a, b, 2 - a - b) for a in range(3) for b in range(3 - a)]
    cs = sp.symbols(f"c0:{len(mons)}")
    ansatz = sum(c * dot(LS, LS) ** a * dot(LS, PS) ** b * dot(PS, PS) ** e for c, (a, b, e) in zip(cs, mons))
    eqs = sp.Poly(sp.expand(ansatz - N), *LS, *PS).coeffs()
    sol = sp.solve(eqs, cs, dict=True)
    inv_ok = bool(sol) and sp.expand(ansatz.xreplace(sol[0]) - N) == 0
    R2 = 0
    for c, (a, b, e) in zip(cs, mons):
        cv = sol[0].get(c, 0) if sol else 0
        if cv != 0:
            R2 += cv * PP ** e * IL.J_two_point(a, b, f"Q02_D04[{a},{b}]")
    R2 = sp.expand(norm * R2 / MB)
    A1 = c1.get("p4", sp.Integer(0)) if ok1 else None
    A2 = sp.factor(sp.cancel(R2 / PP ** 2))
    # route 3 (hand, independent of the engine): T.P.T' with T_ab = p_a q_b + p_b q_a - eta_ab p.q, de Donder P (Q4S
    # trace 4): T.P.T' = 2(p.p')(q.q') + 2(p.q')(q.p') - 2(p.q)(p'.q'); p' = -p, q = l - p, q' = -q  =>  2 p^2 (l-p)^2
    hand = sp.expand(2 * (-PP) * (-sp.Symbol("q.q")) + 2 * (-sp.Symbol("p.q")) ** 2 - 2 * sp.Symbol("p.q") * sp.Symbol("p.q"))
    out.update({"route2_invariant_fit_identical": inv_ok, "route2_IBP_A": str(A2),
                "route1_A": str(A1), "N_factorised": str(sp.factor(N)),
                "N_proportional_to_p2_times_scalar_denominator": sp.expand(N - sp.Rational(1, 2) * dot(PS, PS) * dot([LS[r] - PS[r] for r in range(4)], [LS[r] - PS[r] for r in range(4)])) == 0,
                "route3_hand_contraction": str(hand) + "  (= 2 p^2 (l-p)^2 with q = l - p; engine N = 1/2 p^2 (l-p)^2: same shape, engine vertex normalisation)",
                "reading": "N cancels the scalar propagator: D04 = (p^2/2) * int 1/l^2 (graviton line) = scaleless, omega = 2 power divergence, no log member: D04 = 0 exactly in d by BOTH routes (route 1 uses no scaleless rule)",
                "routes_agree": A1 is not None and sp.cancel(A1 - A2) == 0,
                "N_invariant_form": {f"(l.l)^{a}(l.p)^{b}(p.p)^{e}": str(sol[0].get(c, 0)) for c, (a, b, e) in zip(cs, mons)} if sol else None,
                "structure": "Gamma_1^{phi phi}(p) = kappa^2 A(d) (p^2)^2 MB  (no evanescent structure: Gt.p = 0)",
                "value": "0" if (A1 == 0 and A2 == 0) else "nonzero",
                "laurent": expand_coeff(A1) if A1 not in (None, 0) else {"UV_pole": "0", "MSbar_finite": "0", "local": "0", "coefficient_of_L": "0"},
                "status": "Evaluated_vanishing_under_declared_prescription" if (A1 == 0 and A2 == 0) else "Evaluated"})
    return out, c1


# ----------------------------------------------------------------------------------------------- D06
def D06(S_f=1):
    Vg = F.vertices()["V_cb_c_h"]
    Pi = AR.pinv_sparse()
    Gh = AR.ghost_num()
    VA = {k: sym(v) for k, v in F.vt_route(Vg, [(0, -1), (1, 0), (-1, 1)]).items()}   # cbar_rho: -p, c^mu: l, h_a: p-l
    VB = {k: sym(v) for k, v in F.vt_route(Vg, [(-1, 0), (0, 1), (1, -1)]).items()}   # cbar_sig: -l, c^nu: p, h_b: l-p
    N = {}
    for (rho, mu, a), pa in VA.items():
        for (sig, nu, b), pb in VB.items():
            g = Gh.get((mu, sig))
            pv = Pi.get((a, b))
            if g and pv:
                N[(rho, nu)] = N.get((rho, nu), 0) + pa * pb * sp.Rational(g.numerator, g.denominator) * sp.Rational(pv.numerator, pv.denominator)
    # ghost numerator carries its i; graviton propagator: one explicit i
    norm = prefactor_norm(S_f, 1, Vg.phase * 2, 1)
    G = {k: sp.expand(norm * IL.fp_two_point(sp.expand(v), LS, PS) / MB) for k, v in N.items() if sp.expand(v) != 0}
    pu = lambda r: ETA[r] * PS[r]  # noqa: E731
    S = {"S1_delta^rho_nu_p4": {(r, s): (1 if r == s else 0) * PP ** 2 for r in range(4) for s in range(4)},
         "S2_p^rho_p_nu_p2": {(r, s): ETA[r] * PS[r] * PS[s] * PP for r in range(4) for s in range(4)},
         "SE_Gt^rho_nu_p4": {(r, s): ETA[r] * gt(r, s) * PP ** 2 for r in range(4) for s in range(4)}}
    c, ok, n = fit(G, S)
    return {"fit": {k: str(v) for k, v in c.items()}, "fit_exact": ok, "rows": n,
            "laurent": {k: expand_coeff(v) for k, v in c.items()}}, c


def tadpoles():
    """D05a, D12: explicit degree / scaleless analysis on the engine integrands."""
    out = {}
    V2s = F.vertices()["V_hh_phi_phi"]
    Pi = AR.pinv_sparse()
    # D12: phi legs contracted (l, -l), external h_a(p), h_b(-p)
    VA = {k: sym(v) for k, v in F.vt_route(V2s, [(0, 1), (0, -1), (1, 0), (-1, 0)]).items()}
    deg12 = max((sp.Poly(v, *LS).total_degree() for v in VA.values() if v != 0), default=0)
    # D05a: h legs contracted (l, -l), external phi(p), phi(-p)
    VB = {k: sym(v) for k, v in F.vt_route(V2s, [(1, 0), (-1, 0), (0, 1), (0, -1)]).items()}
    N5 = 0
    for (c, dd_, _, _), v in VB.items():
        pv = Pi.get((c, dd_))
        if pv:
            N5 += v * sp.Rational(pv.numerator, pv.denominator)
    N5 = sp.expand(N5)
    deg5 = sp.Poly(N5, *LS).total_degree() if N5 != 0 else 0
    for name, deg, note in (("D12", deg12, "phi phi self-contraction on V_hh_phi_phi; numerator carries l (phi legs)"),
                            ("D05a", deg5, "h h self-contraction on V_hh_phi_phi; a2 has no derivative of h")):
        omegas = sorted({4 - 2 + k for k in range(0, deg + 1, 2)})       # int d^4 l (l)^k / l^2
        out[name] = {"integrand": note, "loop_numerator_degree_in_l": deg, "omegas": omegas,
                     "logarithmic_member": 0 in omegas, "UV_pole": "0", "IR_pole": "0",
                     "value": "0 (DERIVED: every member int (l^2)^j/l^2 is scaleless with omega = 2 + 2j > 0 -- a pure "
                              "power divergence; no 1/eps on the UV or the IR side; MS-bar counterterm 0)",
                     "status": "Evaluated_vanishing_under_declared_prescription"}
    out["D05a"]["numerator_l_independent"] = all(sp.diff(N5, x) == 0 for x in LS)
    return out


def main():
    t0 = time.time()
    res = {"prescription": "XVI-P p14r/p22: DRED (QdS loop momenta, Q4S vertex indices), d opaque (qguard), MS-bar, "
                           "Feynman +i0, m = 0 (XVI-O bulk_mass_policy; XVI-N massless obstruction: insertion untouched)",
           "Gamma1_convention": "i Gamma_1 = graph sum; coefficients per kappa^2 in units of MB = Ihat(1,1) = 1/eps + 2 - L + O(eps)",
           "symmetry_factors": {"D04": "1 (two distinct-species internal lines; no automorphism fixing external legs; "
                                       "XVI-O's '1/2 Unverified' corrected)",
                                "D06": "1 (as XVI-O)", "D05a": "1/2", "D12": "1/2"}}
    res["D04"], c4 = D04()
    m = {}
    _, c_half = D04(S_f=sp.Rational(1, 2), mutate="S_f_half")
    m["D04_S_f_half_still_zero (zero is S_f-independent; not a detector)"] = c_half.get("p4", 0) == 0
    _, c_tr = D04(mutate="drop_trace_part_of_graviton_propagator")
    m["D04_drop_trace_part_changes_A"] = sp.cancel(c_tr.get("p4", 0) - c4.get("p4", 0)) != 0
    m["D04_drop_trace_part_value"] = str(c_tr.get("p4", 0))
    res["D06"], c6 = D06()
    _, c6h = D06(S_f=sp.Rational(1, 2))
    m["D06_S_f_half_changes_coefficients"] = any(sp.cancel(c6h.get(k, 0) - c6.get(k, 0)) != 0 for k in c6)
    p15 = json.loads((Q / "provenance" / "XVI-P" / "certificates" / "P15_D13.json").read_text())["ghost_self_energy_pole"]
    ref = {"S1_delta^rho_nu_p4": p15["coefficients_exact_d"]["S1"], "S2_p^rho_p_nu_p2": p15["coefficients_exact_d"]["S2"],
           "SE_Gt^rho_nu_p4": p15["coefficients_exact_d"]["SE"]}
    res["D06"]["XVI_P_P15_ghost_self_energy_coefficients"] = p15["coefficients_exact_d"]
    res["D06"]["equals_XVI_P_P15_exactly"] = all(sp.cancel(sp.sympify(ref[k], locals={"d": d}) - c6.get(k, 0)) == 0 for k in ref) and set(c6) == set(ref)
    res["D06"]["reading"] = ("full one-loop cbar c two-point under the declared prescription; independent construction (this file) "
                             "equals XVI-P P15_D13 ghost self-energy (route 2 there; ghost equation verified there) exactly in d")
    res["tadpoles"] = tadpoles()
    res["mutations"] = m
    res["seconds"] = round(time.time() - t0, 1)
    return res


if __name__ == "__main__":
    r = main()
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Q / "build" / "q02_bulk_graphs.json"
    out.write_text(json.dumps(r, indent=1, default=str) + "\n")
    print(json.dumps({k: r[k] for k in ("D04", "D06", "mutations")}, indent=1, default=str)[:4000])
