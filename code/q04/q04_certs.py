#!/usr/bin/env python3
"""XVI-Q Q04 certificates (P.6 entry, one-loop kappa^3; honest partial):
  Q04_P6_inventory.json  independent enumeration (15 loop graphs == XVI-P P23 + 3 kappa^3 tadpoles), power counting,
                         kappa^3 vertices (V_hhh_phi_phi exact: two routes + mutation; V5_grav enters only a scaleless tadpole)
  Q04_P6_identity.json   kappa^3 identity ledger: terms, gh1 three-field rows, inputs per row, CT-domain structural counts
  Q04_P6_eval.json       own UV-pole engine validated against XVI-P D11 (all 14 structures incl. evanescent) + mutation;
                         scalar-loop hhh graphs G14/G15 (UV pole residue, exact); matter Ward identity at three points
                         (pole residue) + mutations; graph table with every other graph explicitly Missing
  Q04_P6_status.json     outcome per stage; defect/CT/closure status"""
from __future__ import annotations
import hashlib, json, sys, time  # noqa: E401
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

Q = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Q / "code" / "q04"))
CERT = Q / "certificates"


def _sub4(e, why):
    import dred
    import sympy as sp
    return sp.expand(dred.substitute_d(e, 4, why))


def task_d11(trace_coef="1/2"):
    import sympy as sp
    import uvpole as U
    U.TRACE_COEF = sp.Rational(trace_coef)
    p, l = U.vec("p1"), U.vec("l")
    e1, e2 = U.rank1("a1"), U.rank1("a2")
    G2 = U.loop_pole(U.V_hphiphi(e1, l, -(l + p)) * U.V_hphiphi(e2, l + p, -l), [U.V(), p], 4, sp.Rational(1, 2))
    mine = _sub4(G2, "Q04 D11 validation: MS-bar residue")
    r = json.loads((Q / "provenance/XVI-P/certificates/P18_D11.json").read_text())["result"]
    H, T = U.H, U.T
    idx = {"m": 1, "n": 1, "r": 2, "s": 2}
    tot = 0
    for k, v in r.items():
        parts = k.split(" ")
        e = int(parts[-1][-1]) if parts[-1].startswith("(p.p)") else 0
        parts = [f for f in parts if not f.startswith("(p.p)")]
        val = H("p1", "p1") ** e
        orbit = set()
        for s1 in (0, 1):
            for s2 in (0, 1):
                mp = {"m": "nm"[s1], "n": "mn"[s1], "r": "sr"[s2], "s": "rs"[s2]}
                ren = []
                for f in parts:
                    nm, ix = f.split("_")
                    nix = "".join(mp[c] for c in ix)
                    ren.append(nm + "_" + ("".join(sorted(nix)) if nm in ("g4", "gt") else nix))
                orbit.add(tuple(sorted(ren)))
        for f in parts:
            nm, ix = f.split("_")
            a, b = (f"a{idx[ix[0]]}", f"a{idx[ix[1]]}") if len(ix) == 2 else (f"a{idx[ix]}", None)
            val *= {"g4": lambda: H(a, b) + T(a, b), "gt": lambda: T(a, b), "p": lambda: H(a, "p1")}[nm]()
        tot += len(orbit) * sp.sympify(v["UV_pole"]) * val
    ratio = sp.cancel(mine / sp.expand(tot))
    return {"trace_coef": trace_coef, "ratio_mine_over_XVI_P": str(ratio), "constant_ratio": not ratio.free_symbols,
            "n_structures": len(r), "gt_reading": "tilde (evanescent Q(4-d)S metric); g4 = hat + tilde; XVI-P 'upper' "
            "structures weighted by their (m<->n, r<->s) orbit size"}


def task_ward(lie_sign=1, drop_bubbles=False):
    import sympy as sp
    import scalar_loop3 as S
    w, g3 = S.ward(drop_bubbles=drop_bubbles, lie_sign=lie_sign)
    w4 = _sub4(w, "Q04 matter Ward identity: pole residue")
    Ts = [s for s in w.free_symbols if s.name.startswith("T_")]
    return {"lie_sign": lie_sign, "drop_bubbles": drop_bubbles,
            "residue_d4_zero": w4 == 0, "residue_d4_4D_part_zero": sp.expand(w4.subs({s: 0 for s in Ts})) == 0,
            "n_terms_residue_d4": 0 if w4 == 0 else len(sp.Add.make_args(w4)),
            "coefficient_exact_in_d_zero": w == 0, "n_terms_generic_d": len(sp.Add.make_args(w)) if w != 0 else 0}


def task_gamma3():
    import sympy as sp
    import scalar_loop3 as S
    import uvpole as U
    p1, p2 = U.vec("p1"), U.vec("p2")
    g3, parts = S.gamma3([(U.rank1("a1"), p1), (U.rank1("a2"), p2), (U.rank1("a3"), -(p1 + p2))])
    g4 = _sub4(g3, "Q04 scalar-loop hhh: MS-bar residue")
    Ts = [s for s in g4.free_symbols if s.name.startswith("T_")]
    terms = sorted(str(x) for x in sp.Add.make_args(g4))
    per = {k: len(sp.Add.make_args(_sub4(v, "Q04 per-graph residue"))) for k, v in parts.items()}
    return {"residue_d4_n_terms": len(terms), "residue_d4_sha256": hashlib.sha256("\n".join(terms).encode()).hexdigest(),
            "residue_d4_n_evanescent_terms": sum(1 for x in sp.Add.make_args(g4) if x.free_symbols & set(Ts)),
            "residue_d4": str(g4), "per_part_n_terms": per,
            "bose_symmetric_a1_a2": _bose(g4)}


def _bose(g4):
    """swap (a1, p1) <-> (a2, p2) (leg 3 fixed, p3 = -p1-p2 invariant), re-canonicalising H/T symbol names."""
    import sympy as sp
    sw = {"a1": "a2", "a2": "a1", "p1": "p2", "p2": "p1"}
    rep = {}
    for sym in g4.free_symbols:
        k, x, y = sym.name.split("_")
        x, y = sorted((sw.get(x, x), sw.get(y, y)))
        rep[sym] = sp.Symbol(f"{k}_{x}_{y}")
    return sp.expand(g4 - g4.xreplace(rep)) == 0


def main():
    t0 = time.time()
    import p6_inventory as P
    inv = P.inventory()
    inv["V_hhh_phi_phi"] = {"formula": P.X3_formula.__doc__.strip(), "two_routes_agree": P.vertex_hhhphiphi_check(),
                            "mutation_tr_h_cubed_detected": not P.vertex_hhhphiphi_check(mutate=True),
                            "role": "enters only the scaleless tadpoles (E3 tadpoles 2, 3)"}
    inv["V5_grav"] = {"status": "NOT constructed", "why_not_needed": "enters only the self-contracted tadpole (hhh), which "
                      "is scaleless for ANY polynomial vertex; no other one-loop E = 3 graph contains it"}
    (CERT / "Q04_P6_inventory.json").write_text(json.dumps(inv, indent=1) + "\n")
    ident = P.identity_rows()
    ident["counterterm_domain_kappa3"] = P.ct_domain_dims()
    (CERT / "Q04_P6_identity.json").write_text(json.dumps(ident, indent=1) + "\n")
    with ProcessPoolExecutor(max_workers=6) as ex:
        f = {"d11": ex.submit(task_d11), "d11_mut": ex.submit(task_d11, "1"), "ward": ex.submit(task_ward),
             "ward_mut_lie": ex.submit(task_ward, -1), "ward_mut_bub": ex.submit(task_ward, 1, True),
             "g3": ex.submit(task_gamma3)}
        R = {k: v.result() for k, v in f.items()}
    p23 = json.loads((Q / "provenance/XVI-P/certificates/P23_O1_scope.json").read_text())["P6_entry_kappa3_three_point"]["graphs"]
    table = []
    for i, g in enumerate(p23, 1):
        evald = g["vertices"] in (["V_h_phi_phi"] * 3, ["V_h_phi_phi", "V_hh_phi_phi"]) and g["external"] == "h,h,h"
        table.append({"id": f"P6-G{i:02d}", "vertices": g["vertices"], "lines": g["loop_lines"], "external": g["external"],
                      "status": "UV pole evaluated (exact residue, own rules, validated vs D11; Ward-checked); finite part Missing"
                      if evald else "Missing (not evaluated in Q04)"})
    for j, tp in enumerate(inv["E3_tadpoles"], 1):
        table.append({"id": f"P6-T{j}", "vertices": tp["vertices"], "lines": tp["lines"], "external": ",".join(tp["external"]),
                      "status": "derived 0 (scaleless, no log member)"})
    ev = {"engine": "code/q04/uvpole.py (own Feynman rules, symbolic d via qguard, DRED hat/tilde split, IR-rearranged "
                    "Taylor UV pole; units U = i/(16 pi^2 eps); 'i Gamma' convention)",
          "validation_D11": R["d11"], "validation_mutation_trace_coef_1": R["d11_mut"],
          "scalar_loop_hhh": R["g3"],
          "matter_ward_three_point": {"identity": "Gamma3(d0h, e1, e2) + Gamma2(d1h[xi,e1], e2) + Gamma2(d1h[xi,e2], e1) = 0 "
                                                  "(diffeomorphism invariance of the scalar-loop W[g]); row c.h.h, matter subsector",
                                      "declared": R["ward"], "mutation_wrong_Lie_sign": R["ward_mut_lie"],
                                      "mutation_drop_X2_bubbles": R["ward_mut_bub"],
                                      "reading": "the MS-bar pole residue satisfies the identity exactly, incl. the evanescent "
                                                 "(tilde) structures; the d-dependent coefficient differs by O(d-4) x U, i.e. "
                                                 "only finite (non-pole) terms, which Q04 does not evaluate"},
          "graph_table": table, "n_evaluated_pole": sum(1 for x in table if x["status"].startswith("UV pole")),
          "n_missing": sum(1 for x in table if x["status"].startswith("Missing")),
          "n_derived_zero": sum(1 for x in table if x["status"].startswith("derived 0"))}
    (CERT / "Q04_P6_eval.json").write_text(json.dumps(ev, indent=1) + "\n")
    ok_val = R["d11"]["constant_ratio"] and R["d11"]["ratio_mine_over_XVI_P"] == "I" and not R["d11_mut"]["constant_ratio"]
    ok_ward = R["ward"]["residue_d4_zero"] and not R["ward_mut_lie"]["residue_d4_zero"] and not R["ward_mut_bub"]["residue_d4_zero"]
    st = {"landmark": "XVI-Q", "release": "Q04", "P6_outcome": "Partial",
          "stages": {"inventory": "Constructed" if inv["equals_XVI_P_P23_15"] and inv["n_E3_tadpoles"] == 3 else "FAIL",
                     "kappa3_vertices": "Constructed (V_hhh_phi_phi); V5_grav not needed (scaleless only)",
                     "power_counting_kappa3_1pt_2pt_empty": inv["power_counting"]["verified_on_all_enumerated_graphs"],
                     "identity_ledger_rows": f"Constructed ({ident['n_rows']} gh1 three-field rows; "
                                             f"{sum(1 for r in ident['rows'] if r['class'] == 'needs_inputs')} need P.6 inputs)",
                     "counterterm_domain_kappa3": "Partial (structural counts; BRST/sigma_0 reduction not done)",
                     "engine_validated_vs_XVI_P_D11": ok_val,
                     "evaluation": f"Partial: {ev['n_evaluated_pole']} graphs at pole level, {ev['n_derived_zero']} tadpoles "
                                   f"derived 0, {ev['n_missing']} Missing",
                     "matter_subsector_ward_c_h_h_pole": "verified" if ok_ward else "FAIL",
                     "defect_assembly": "Missing (inputs Missing for every needs_inputs row)",
                     "counterterm_solve_kappa3": "Missing", "sigma0_closure_kappa3": "Missing"},
          "rows_status": {r["row"]: ("zero (" + r["why"][:40] + ")" if r["class"] == "zero" else
                                     "aux" if r["class"] == "aux" else "Missing inputs") for r in ident["rows"]},
          "O1_closure_unchanged": "CLOSED (scoped) -- kappa^2 two-point; P.6 is outside that scope",
          "bulk_symmetry_restored": False, "seconds": round(time.time() - t0, 1)}
    (CERT / "Q04_P6_status.json").write_text(json.dumps(st, indent=1) + "\n")
    import guard_scan
    (CERT / "Q04_guard_scan.json").write_text(json.dumps(guard_scan.report(), indent=1) + "\n")
    print(json.dumps({k: v for k, v in st["stages"].items()}, indent=0))


if __name__ == "__main__":
    main()
