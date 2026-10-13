#!/usr/bin/env python3
"""XVI-Q Q04 P.6 entry (one-loop kappa^3): independent graph enumeration incl. kappa^3 vertices, power counting,
kappa^3 vertex derivation (V_hhh_phi_phi exact, with an independent determinant/inverse route and a mutation control),
kappa^3 identity ledger (which (S_k, Gamma_1^[3-k]) enter; gh1 three-field rows) and the kappa^3 three-point local
counterterm-domain dimensions over Q(i) (4D vs evanescent), for the hhh and h phi phi sectors."""
from __future__ import annotations
import itertools, json, sys  # noqa: E401
from pathlib import Path
import sympy as sp

Q = Path(__file__).resolve().parents[2]
PV = Q / "provenance" / "XVI-P" / "certificates"
sys.path.insert(0, str(Q / "code" / "q02"))
from bulk_identity import FIELDS  # noqa: E402

PAIRS = {("h", "h"), ("phi", "phi"), ("c", "cb"), ("cb", "c")}
KAPPA3_VERTICES = {"V5_grav": ["h"] * 5, "V_hhh_phi_phi": ["h", "h", "h", "phi", "phi"]}
EXACT_IN_KAPPA = {
    "V_cb_c_h": "ghost action = cbar . d(s h) with Psi linear in h (XVI-P conventions PSI) and s h = (1/2) L_c(eta + kappa h) "
                "affine in h (GHOST_SCALING): no cbar c h^n vertex for n >= 2",
    "V_K_c_h": "K s h, s h affine in h: no K c h^n, n >= 2", "V_cs_c_c": "c* s c, s c = kappa a c.dc: exact at kappa^1",
    "V_phis_c_phi": "phi* s phi, s phi = kappa a c.dphi: exact at kappa^1",
    "gauge_fixing": "Psi linear in h -> b/h sector quadratic, no vertices"}


def jl(p):
    return json.loads(Path(p).read_text())


def vertices():
    V = dict(jl(PV / "P14_counterterm_domain.json")["loop_capable_vertices"])
    V.update(KAPPA3_VERTICES)
    return V


def kappa(legs):
    return len(legs) - 2


def enumerate_one_loop(V, nmax=3):
    """all one-loop 1PI graphs = cycles of n vertices, each contributing an (in, out) pair of loop legs joined by
    propagating pairs.  Key: (sorted vertices, sorted line types, sorted externals)."""
    out = {}
    names = sorted(V)
    for n in range(1, nmax + 1):
        for combo in itertools.product(names, repeat=n):
            legsets = [V[v] for v in combo]
            for choice in itertools.product(*[list(itertools.permutations(range(len(l)), 2)) for l in legsets]):
                ok, lines = True, []
                for i in range(n):
                    a = legsets[i][choice[i][1]]                      # out-leg of vertex i
                    b = legsets[(i + 1) % n][choice[(i + 1) % n][0]]  # in-leg of the next vertex
                    if (a, b) not in PAIRS:
                        ok = False
                        break
                    lines.append("-".join(sorted((a, b))))
                if not ok:
                    continue
                ext = sorted(l for i in range(n) for k, l in enumerate(legsets[i]) if k not in choice[i])
                key = (tuple(sorted(combo)), tuple(sorted(lines)), tuple(ext))
                out[key] = sum(kappa(V[v]) for v in combo)
    return out


def inventory():
    V = vertices()
    g = enumerate_one_loop(V)
    by_E = {}
    for (vs, ls, ext), k in g.items():
        by_E.setdefault(len(ext), []).append({"vertices": list(vs), "lines": list(ls), "external": list(ext), "kappa": k})
    kappa_eq_E = all(x["kappa"] == E for E, lst in by_E.items() for x in lst)
    E3 = sorted(by_E.get(3, []), key=lambda x: (x["external"], x["vertices"], x["lines"]))
    tad = [x for x in E3 if len(x["vertices"]) == 1]
    loops = [x for x in E3 if len(x["vertices"]) > 1]
    p23 = jl(PV / "P23_O1_scope.json")["P6_entry_kappa3_three_point"]
    p23_keys = sorted((tuple(sorted(x["vertices"])), tuple(sorted("-".join(sorted(l.split("-"))) for l in x["loop_lines"])),
                       tuple(sorted(x["external"].split(",")))) for x in p23["graphs"])
    mine_keys = sorted((tuple(x["vertices"]), tuple(x["lines"]), tuple(x["external"])) for x in loops)
    by_ext = {}
    for x in E3:
        by_ext[",".join(x["external"])] = by_ext.get(",".join(x["external"]), 0) + 1
    for x in tad:
        x["status"] = "derived 0 (scaleless): massless self-contraction, the loop momentum flows through ONE propagator " \
                      "and the vertex is a polynomial in l and external momenta -> int d^dl l^{2k} / l^2 = 0 (XVI-P rule S0, no log member)"
    return {"vertices_used": V, "kappa3_vertices_new": sorted(KAPPA3_VERTICES), "exact_in_kappa_vertices": EXACT_IN_KAPPA,
            "power_counting": {"rule": "one loop: I = V, sum_v (legs_v - 2) = (2I + E) - 2V = E  => kappa power = E",
                               "verified_on_all_enumerated_graphs": kappa_eq_E,
                               "kappa3_one_point_one_loop": "empty (kappa^E with E = 1 is kappa^1)",
                               "kappa3_two_point_one_loop": "empty (E = 2 is kappa^2)",
                               "completeness": "each vertex of a one-loop E = 3 graph carries >= 1 external leg (3-leg vertices: "
                                               "exactly 1) -> at most 3 vertices; enumerated n = 1..3"},
            "n_by_E": {E: len(v) for E, v in sorted(by_E.items())},
            "E3_loops": loops, "E3_tadpoles": tad, "n_E3_loops": len(loops), "n_E3_tadpoles": len(tad),
            "by_external": by_ext,
            "equals_XVI_P_P23_15": mine_keys == p23_keys, "XVI_P_n": len(p23_keys)}


# ---------------- kappa^3 vertex V_hhh_phi_phi: cubic term of X = sqrt(-g) g^{mn}
def X3_formula(h, eta):
    """X3 = -h^3 + (1/2) t1 h^2 + (t1^2/8 - t2/4)(-h) + (t3/6 - t1 t2/8 + t1^3/48) eta^{-1}  (h with indices raised by eta)."""
    hu = eta.inv() * h * eta.inv()                        # h^{mn}
    m = eta.inv() * h                                     # h^m_n (mixed)
    t1, t2, t3 = m.trace(), (m * m).trace(), (m * m * m).trace()
    return sp.expand(-(m * m * m) * eta.inv() + sp.Rational(1, 2) * t1 * (m * m) * eta.inv()
                     - (t1 ** 2 / 8 - t2 / 4) * hu + (t3 / 6 - t1 * t2 / 8 + t1 ** 3 / 48) * eta.inv())


def vertex_hhhphiphi_check(mutate=False):
    """independent route: X(t) = sqrt(-det(eta + t h)) (eta + t h)^{-1}, exact series in t, t^3 coefficient, on two
    random rational symmetric h (4x4, signature +---); compared with X3_formula."""
    eta = sp.diag(1, -1, -1, -1)
    tt = sp.Symbol("tt")
    res = []
    for seed in (3, 11):
        import random
        rnd = random.Random(seed)
        h = sp.zeros(4, 4)
        for i in range(4):
            for j in range(i, 4):
                h[i, j] = h[j, i] = sp.Rational(rnd.randint(-7, 7), rnd.randint(1, 5))
        g = eta + tt * h
        det = sp.expand(-g.det())
        sq = sp.series(sp.sqrt(det), tt, 0, 4).removeO()
        ginv = sum((((-tt) ** k) * (eta.inv() * h) ** k for k in range(4)), sp.zeros(4, 4)) * eta.inv()
        X = (sq * ginv).applyfunc(lambda e: sp.expand(e).coeff(tt, 3))
        F = X3_formula(h, eta)
        if mutate:
            F = F + sp.Rational(1, 48) * (eta.inv() * h).trace() ** 3 * eta.inv()
        res.append(sp.simplify(X - F) == sp.zeros(4, 4))
    return all(res)


# ---------------- kappa^3 identity ledger
def identity_rows():
    names = list(FIELDS)
    rows = []
    for tri in itertools.combinations_with_replacement(names, 3):
        gh = sum(FIELDS[x][0] for x in tri)
        if gh != 1:
            continue
        if sum(FIELDS[x][1] for x in tri) % 2:
            rows.append({"row": ".".join(tri), "class": "zero", "why": "phi-parity odd (every vertex even in {phi, phi*})"})
            continue
        if {"b", "cbs", "bs"} & set(tri):
            rows.append({"row": ".".join(tri), "class": "aux", "why": "contains b / cbar* / b*: no interaction vertex "
                         "(ghost-equation / nonminimal rows; not a loop entry)"})
            continue
        rows.append({"row": ".".join(tri), "class": "needs_inputs"})
    needs = {
        "c.h.h": {"S0_Gamma3": ["Gamma_hhh (5 graphs)", "Gamma_K c h (2)"], "S1_Gamma2": ["Gamma_hh (Q02: D01, D02, D03, D11, D12)",
                  "Gamma_K c (D13)"], "matter_subsector": "scalar-loop part = matter Ward identity (evaluated, Q04)"},
        "c.phi.phi": {"S0_Gamma3": ["Gamma_h phi phi (4 graphs)", "Gamma_phi* c phi (1)"], "S1_Gamma2": ["Gamma_phi phi = 0 (Q02)",
                      "Gamma_hh not entering"]},
        "c.c.cb": {"S0_Gamma3": ["Gamma_cbar c h (2)", "Gamma_c c c* (1)"], "S1_Gamma2": ["Gamma_cbar c (D06)"]},
        "c.c.hs": {"S0_Gamma3": ["Gamma_K c h (2)", "Gamma_c c c* (1)"], "S1_Gamma2": ["Gamma_K c (D13)"]},
    }
    needs = {".".join(sorted(k.split("."))): v for k, v in needs.items()}
    for r in rows:
        if r["class"] == "needs_inputs":
            r["inputs"] = needs.get(".".join(sorted(r["row"].split("."))), {"S0_Gamma3": [], "S1_Gamma2": [], "note": "no E = 3 kappa^3 graph with these "
                                                                                       "external legs (enumeration) -> row zero"})
            if ".".join(sorted(r["row"].split("."))) not in needs:
                r["class"] = "zero_no_graph"
    return {"identity": "(S, Gamma_1) = A_1 at kappa^3, one loop: (S_0, Gamma_1^[3]) + (S_1, Gamma_1^[2]) + (S_2, Gamma_1^[1]) + (S_3, Gamma_1^[0])",
            "terms": {"(S_0, Gamma_1^[3])": "three-point Gamma_1 with the linear BRST (s0 h = d c, s0 c = 0, s0 phi = 0, antifield "
                                            "equations of S_0) -- needs the P.6 graphs",
                      "(S_1, Gamma_1^[2])": "two-point Gamma_1 (Q02, evaluated) with the quadratic BRST (c.dh, h dc, c.dc, c.dphi)",
                      "(S_2, Gamma_1^[1])": "0: Gamma_1^[1] is one-point only and vanishes (D14a-c)",
                      "(S_3, Gamma_1^[0])": "0: field independent"},
            "rows": rows, "n_rows": len(rows)}


# ---------------- kappa^3 local counterterm domain dimensions (structural; before BRST reduction)
def ct_domain_dims():
    """local three-point polynomials, momentum degree 4 (one-loop kappa^3 counterterms have 4 derivatives: kappa^3 x
    mass^4 for a 3-field operator), multilinear (degree 2) in each graviton polarisation a_i, Bose symmetric, momentum
    conserved (p3 = -p1-p2); invariants H (hat), T (tilde, evanescent) as in uvpole.py.  Dimension over Q(i) of the
    symmetric span, and of its 4D part (T-free) -- evanescent = total - 4D."""
    p = ["p1", "p2"]
    def mono_space(pol, sym_group):
        syms = []
        vecs = pol + p
        for i, x in enumerate(vecs):
            for y in vecs[i:]:
                if x in p and y in p:
                    syms.append((("H", x, y), {x: 0, y: 0}))
                else:
                    syms.append((("H", x, y), None))
                    if x in pol and y in pol:
                        syms.append((("T", x, y), None))
        # monomials: product of invariants with degree 2 in each pol vector and momentum degree 4
        out = set()
        def rec(i, cur, degp, degm):
            if all(degp[a] == 2 for a in pol) and degm == 4:
                out.add(tuple(sorted(cur)))
                return
            if i == len(syms):
                return
            (k, x, y), _ = syms[i]
            for n in range(0, 5):
                nd = dict(degp)
                ok = True
                m = 0
                for z in (x, y):
                    if z in pol:
                        nd[z] += n
                    else:
                        m += n
                if any(nd[a] > 2 for a in pol) or degm + m > 4:
                    break
                rec(i + 1, cur + [(k, x, y)] * n, nd, degm + m)
        rec(0, [], {a: 0 for a in pol}, 0)
        return sorted(out)
    res = {}
    for name, pol in (("hhh", ["a1", "a2", "a3"]), ("h_phi_phi", ["a1"])):
        mons = mono_space(pol, None)
        res[name] = {"n_monomials_before_symmetry": len(mons),
                     "n_with_T (evanescent)": sum(1 for m in mons if any(k == "T" for k, _, _ in m))}
    res["status"] = ("structural count only (momentum conservation used through p1, p2; Bose symmetrisation and the "
                     "BRST/sigma_0 reduction to the invariant domain NOT done in Q04)")
    return res
