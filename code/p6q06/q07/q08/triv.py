#!/usr/bin/env python3
"""XVI-Q Q08 goals 1-3: structure of the local evanescent kappa^3 counterterm X (three-point rows).
 A x = -b : rows (c.c.K, c.c.cbar, h.h.c) linear in the X coordinates (Q07 ansatz incl. Bose symmetry of X_hhh).
 (1) ghost equation: add X_{cbar ch} = rho X_{Kch}[K -> T(cbar)], rho = -2 (fixed by the loop poles, geq.py).
 (3) homogeneous space N = ker A; BRST-trivial directions sigma_lin(Y) for gh -1 local Y at three points:
     F: Y = K.F(h,h)        -> X_hhh = sum_i dS2/dh_i . F_i,      X_Khc ~ K.F(s0h c, h)
     G: Y = c*.G(c,h)       -> X_cbar ch ~ cbar.D[s0h G],          X_Khc ~ K.s0h(G)
     P: Y = cbar.Phi(h,h)   -> X_hhh ~ sigma_F F(h_i).Phi_i,       X_cbar ch ~ cbar.Phi(s0h c, h)
     the relative normalisation lambda of the two pieces of each family is not derived but FITTED (one constant per
     family, required to be the same for every monomial: a nontrivial check); the triviality claim is conditional
     on that.  Output: dims of N, of N with the ghost equation, of the trivial span inside N, remainder."""
from __future__ import annotations
import json, sys, time  # noqa: E401
from pathlib import Path
import sympy as sp
from sympy.polys.matrices import DomainMatrix

HERE = Path(__file__).resolve().parent
Q = HERE.parents[3]
sys.path.insert(0, str(HERE.parents[2] / "p6q05"))
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import rows as R  # noqa: E402
import gravity as G  # noqa: E402
import uvpole as U  # noqa: E402
import hhc  # noqa: E402
import evfit2 as E2  # noqa: E402
from evfit import ansatz, monomials  # noqa: E402
from geq import Tmap  # noqa: E402
LOG = Q / "build" / "logs"
I = sp.I
RHO = -2


def coords(expr, ms, tag):
    idx = {m: i for i, m in enumerate(ms)}
    v = [0] * len(ms)
    out = 0
    for t in sp.Add.make_args(sp.expand(expr)):
        if t == 0:
            continue
        c, m = t.as_independent(*[s for s in t.free_symbols if s.name[:2] in ("H_", "T_")], as_Add=False)
        if m in idx:
            v[idx[m]] += c / I
        else:
            out += t
    return v, out


def build():
    XK, xk = ansatz({"aK": 2, "gc1": 1, "a3": 2}, 3, "xk")
    XC, xc = ansatz({"gb2": 1, "gc1": 1, "a3": 2}, 4, "xc")
    XS, xs = ansatz({"z": 1, "gc1": 1, "gc2": 1}, 3, "xs")
    XH, xh = ansatz({"a1": 2, "a2": 2, "a3": 2}, 4, "xh")
    sig = sp.sympify(json.loads((LOG / "q07_hhc_tree.json").read_text())["sigma_F"])
    rows = [E2.ghost_row_total("ccK", XK, XS), E2.ghost_row_total("ccb", XC, XS), E2.hhc_ansatz_total(XH, XK, XC, sig)]
    xall = xk + xc + xs + xh
    X = set(xall)
    eqs = []
    for r in rows:
        gens = sorted([s for s in r.free_symbols if s not in X], key=str)
        eqs += sp.Poly(r, *gens).coeffs()
    beq = []
    for mp in ({"p1": R.P2, "p2": R.P1, "a1": "a2", "a2": "a1"}, {"p1": R.P3, "p2": R.P2, "a1": "a3", "a3": "a1"}):
        d = sp.expand(XH - R.remap(XH, mp))
        gens = sorted([s for s in d.free_symbols if s not in X], key=str)
        beq += sp.Poly(d, *gens).coeffs()
    geq = sp.expand(XC - RHO * Tmap(XK))
    gens = sorted([s for s in geq.free_symbols if s not in X], key=str)
    geqs = sp.Poly(geq, *gens).coeffs()
    A, b = sp.linear_eq_to_matrix(eqs + beq, xall)
    Ag, bg = sp.linear_eq_to_matrix(geqs, xall)
    ms = {"K": monomials({"aK": 2, "gc1": 1, "a3": 2}, 3), "C": monomials({"gb2": 1, "gc1": 1, "a3": 2}, 4),
          "S": monomials({"z": 1, "gc1": 1, "gc2": 1}, 3), "H": monomials({"a1": 2, "a2": 2, "a3": 2}, 4)}
    return A, b, Ag, bg, ms, sig


def fam_vectors(ms, sig):
    """for each family and monomial: (piece1 coords, piece2 coords) as full-length vectors."""
    nK, nC, nS, nH = (len(ms[k]) for k in "KCSH")
    off = {"K": 0, "C": nK, "S": nK + nC, "H": nK + nC + nS}
    N = off["H"] + nH

    def full(sector, v):
        w = [0] * N
        for i, c in enumerate(v):
            w[off[sector] + i] = c
        return w

    fams = {}
    leak = {}
    # F family: F(K at p1; a2 at p2, a3 at p3), symmetrised in the two h legs
    Fm = monomials({"aK": 2, "a2": 2, "a3": 2}, 2)
    sw23 = {"p1": R.P1, "p2": R.P3, "a2": "a3", "a3": "a2"}
    rel = {"1": {}, "2": {"p1": R.P2, "p2": R.P1, "a2": "a1"}, "3": {"p1": R.P3, "p2": R.P1, "a2": "a1", "a3": "a2"}}
    mom = {"1": hhc.p1, "2": hhc.p2, "3": hhc.p3}
    out = []
    lk = 0
    for m in Fm:
        f = sp.expand(m + R.remap(m, sw23))
        if f == 0:
            continue
        xh_ = 0
        for i in "123":
            fi = R.remap(f, rel[i]) if rel[i] else f
            fi = R.polarise(fi, "aK", "KU", "KW")
            eh, gf = G.S2_EH_gf(hhc.ext(f"a{i}"), G.probe_mat("E"), mom[i])
            xh_ += R.strip_t(R.contract_id(fi, I * (eh + gf), "K", "E", "h"))
        fk = R.remap(f, {"p1": R.P2, "p2": R.P1})                 # K at p2, a2-leg at p1, a3 at p3
        xk_ = R.strip_t(-I * R.polarise(fk, "a2", "p1", "gc1"))
        vh, o1 = coords(xh_, ms["H"], "H")
        vk, o2 = coords(xk_, ms["K"], "K")
        lk += (o1 != 0) + (o2 != 0)
        out.append((full("H", vh), full("K", vk)))
    fams["F"] = out
    leak["F"] = lk
    # G family: G(gv; gc1 at p1, a3 at p3), vector index gv, momentum degree 2
    out = []
    lk = 0
    for m in monomials({"gv": 1, "gc1": 1, "a3": 2}, 2):
        xc_ = -sp.Rational(1, 2) * U.H("p2", "p2") * hhc.subst_vec(m, "gv", {"gb2": 1})
        xk_ = -I * U.H("aK", "p2") * hhc.subst_vec(m, "gv", {"aK": 1})
        vc, o1 = coords(sp.expand(xc_), ms["C"], "C")
        vk, o2 = coords(sp.expand(xk_), ms["K"], "K")
        lk += (o1 != 0) + (o2 != 0)
        out.append((full("C", vc), full("K", vk)))
    fams["G"] = out
    leak["G"] = lk
    # P family: Phi(gv; a2 at p2, a3 at p3), symmetrised, momentum degree 3
    out = []
    lk = 0
    relv = {"1": {}, "2": {"p1": R.P2, "p2": R.P1, "a2": "a1"}, "3": {"p1": R.P3, "p2": R.P1, "a2": "a1", "a3": "a2"}}
    for m in monomials({"gv": 1, "a2": 2, "a3": 2}, 3):
        f = sp.expand(m + R.remap(m, sw23))
        if f == 0:
            continue
        xh_ = 0
        for i in "123":
            fi = R.remap(f, relv[i]) if relv[i] else f
            D = G.Dvec(hhc.ext(f"a{i}"), U.V({"p1": 1}) if i == "1" else (hhc.F2 if i == "2" else hhc.F3))
            xh_ += sig * hhc.subst_vec(fi, "gv", D)
        fc = R.remap(f, {"p1": R.P2, "p2": R.P1})                # cbar at p2, a2-leg at p1, a3 at p3
        fc = hhc.subst_vec(fc, "gv", {"gb2": 1})
        xc_ = R.strip_t(-I * R.polarise(sp.expand(fc), "a2", "p1", "gc1"))
        vh, o1 = coords(sp.expand(R.strip_t(xh_)), ms["H"], "H")
        vc, o2 = coords(xc_, ms["C"], "C")
        lk += (o1 != 0) + (o2 != 0)
        out.append((full("H", vh), full("C", vc)))
    fams["P"] = out
    leak["P"] = lk
    return fams, leak


def main():
    t0 = time.time()
    A, b, Ag, bg, ms, sig = build()
    dA = DomainMatrix.from_Matrix(A)
    full = DomainMatrix.from_Matrix(A.col_join(Ag))
    n = A.shape[1]
    rA = dA.rank()
    rAg = full.rank()
    res = {"n_unknowns": n, "rank_rows_bose": rA, "dim_N": n - rA, "rank_with_ghost_eq": rAg,
           "dim_N_with_ghost_eq": n - rAg}
    # (1) solvability with the ghost equation
    aug = DomainMatrix.from_Matrix(A.col_join(Ag).row_join(b.col_join(bg)))
    res["solvable_with_ghost_equation"] = aug.rank() == rAg
    # (3) trivial families
    fams, leak = fam_vectors(ms, sig)
    res["family_leak_terms_outside_ansatz"] = leak
    import pickle
    (LOG / "q08_triv.pkl").write_bytes(pickle.dumps({"A": A, "Ag": Ag, "b": b, "bg": bg, "fams": fams}))
    cands = [0, 1, -1, 2, -2, sp.Rational(1, 2), -sp.Rational(1, 2), I, -I, 2 * I, -2 * I, I / 2, -I / 2, 4, -4]
    lam = {}
    best = {}
    for name, pairs in fams.items():
        P = sp.Matrix.hstack(*[sp.Matrix(p) for p, q in pairs])
        Qm = sp.Matrix.hstack(*[sp.Matrix(q) for p, q in pairs])
        S = P.row_join(Qm)
        rS = DomainMatrix.from_Matrix(S).rank()
        rAS = DomainMatrix.from_Matrix(A * S).rank()
        sc = {}
        for l in cands:
            M = P + l * Qm
            sc[str(l)] = DomainMatrix.from_Matrix(M).rank() - DomainMatrix.from_Matrix(A * M).rank()
        lb = max(sc, key=lambda k: sc[k])
        lam[name] = {"n_generators": len(pairs), "dim_span_both_pieces": rS,
                     "dim_N_cap_span_both_pieces_upper_bound": rS - rAS,
                     "best_uniform_lambda": lb, "dim_N_cap_family_at_best_lambda": sc[lb], "scan": sc}
        best[name] = P + sp.sympify(lb) * Qm
    res["families"] = lam
    Ball = sp.Matrix.hstack(*best.values())
    Sall = sp.Matrix.hstack(*[sp.Matrix(v) for pr in fams.values() for pq in pr for v in pq])
    res["dim_N_cap_all_families_upper_bound"] = DomainMatrix.from_Matrix(Sall).rank() - DomainMatrix.from_Matrix(A * Sall).rank()
    res["dim_N_cap_all_families_best_lambda"] = DomainMatrix.from_Matrix(Ball).rank() - DomainMatrix.from_Matrix(A * Ball).rank()
    res["dim_N_mod_trivial_lower"] = res["dim_N"] - res["dim_N_cap_all_families_upper_bound"]
    res["seconds"] = round(time.time() - t0, 1)
    (LOG / "q08_triv.json").write_text(json.dumps(res, indent=1) + "\n")
    print(res)


if __name__ == "__main__":
    main()
