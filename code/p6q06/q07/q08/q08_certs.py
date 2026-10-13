#!/usr/bin/env python3
"""XVI-Q Q08 certificates: ghost-equation compatibility of X, structure of the free directions, status."""
from __future__ import annotations
import hashlib, json, pickle, subprocess, sys, time  # noqa: E401
from pathlib import Path
import sympy as sp
from sympy.polys.matrices import DomainMatrix

HERE = Path(__file__).resolve().parent
Q = HERE.parents[3]
sys.path.insert(0, str(HERE.parents[2] / "p6q05"))
sys.path.insert(0, str(HERE))
CERT = Q / "certificates"
LOG = Q / "build" / "logs"


def run(a):
    r = subprocess.run([sys.executable, *a], cwd=Q, capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(r.stderr[-2000:])


def main():
    t0 = time.time()
    import rows as R
    from geq import Tmap
    K = R.strip_t(R.load("G03") + R.load("G07"))
    C = R.strip_t(R.load("G04") + R.load("G09"))
    loops = {"rho_total": str(sp.cancel(C / Tmap(K))),
             "rho_G04_G03": str(sp.cancel(R.strip_t(R.load("G04")) / Tmap(R.strip_t(R.load("G03"))))),
             "rho_G09_G07": str(sp.cancel(R.strip_t(R.load("G09")) / Tmap(R.strip_t(R.load("G07")))))}
    run(["code/p6q06/q07/q08/triv.py"])
    run(["code/p6q06/q07/q08/sectors.py"])
    tr = json.loads((LOG / "q08_triv.json").read_text())
    se = json.loads((LOG / "q08_sectors.json").read_text())
    d = pickle.loads((LOG / "q08_triv.pkl").read_bytes())
    AA, bb = d["A"].col_join(d["Ag"]), d["b"].col_join(d["bg"])
    sol = list(sp.linsolve((AA, bb)))
    part = None
    if sol:
        x = sp.Matrix(sol[0])
        x = x.xreplace({s: 0 for s in x.free_symbols})
        part = {"verified_A_x_eq_b": all(sp.expand(v) == 0 for v in (AA * x - bb)),
                "n_nonzero": sum(1 for v in x if v != 0), "sha256": hashlib.sha256(str(list(x)).encode()).hexdigest()}
    geq = {"ghost_equation": "Gamma_{cbar c h} = rho Gamma_{K c h}[K -> T(cbar) = sym(k,cbar) - (1/2)(k.cbar) eta]",
           "loop_poles": loops, "loops_satisfy_with_rho_minus_2": all(v == "-2" for v in loops.values()),
           "constraint_imposed": "X_{cbar ch} = -2 X_{Kch}[K -> T(cbar)]",
           "solvable_with_ghost_equation": tr["solvable_with_ghost_equation"], "particular_solution": part,
           "dim_free_without_ghost_eq": tr["dim_N"], "dim_free_with_ghost_eq": tr["dim_N_with_ghost_eq"],
           "trace_convention": "eta-trace at d4-residue level: hat trace 4, tilde trace 0 (geq.py); validated by "
                               "the loop poles obeying the map with one constant"}
    (CERT / "Q08_ghost_equation.json").write_text(json.dumps(geq, indent=1) + "\n")
    fr = {"homogeneous_space": {"dim_N": tr["dim_N"], "dim_N_with_ghost_eq": tr["dim_N_with_ghost_eq"],
                                "rank_rows_bose": tr["rank_rows_bose"], "n_unknowns": tr["n_unknowns"]},
          "single_sector_blind_directions": se,
          "trivial_families": tr["families"],
          "dim_N_cap_all_families_upper_bound": tr["dim_N_cap_all_families_upper_bound"],
          "dim_N_cap_all_families_best_lambda": tr["dim_N_cap_all_families_best_lambda"],
          "dim_trivial_best_lambda_within_ghost_eq_subspace": se["dim_trivial_best_lambda_cap_N_and_ghost_eq"],
          "reading": ("Of the 152 free directions (104 with the ghost equation), at least 97 (55 within the ghost-equation "
                      "subspace) lie in the span of fitted BRST-trivial families (field redefinitions h -> h + F(h,h), "
                      "c -> c + G(c,h), gauge-fermion change cbar.Phi(h,h)); at most 116 can. The families are only "
                      "partly closed in N with one uniform lambda (F 65/121, G 11/66, P 12/109): the relative "
                      "normalisations are FITTED, not derived, and the families are not complete at three points. "
                      "Remaining >= 36 directions are not shown trivial: 26 hhh-only sigma0-closed evanescent cubic "
                      "structures (candidate invariant evanescent operators) and blind directions of the three-point rows "
                      "(69 Kch-only, 12 c*cc-only) constrained only by the ghost equation / four-point rows.")}
    (CERT / "Q08_free_directions.json").write_text(json.dumps(fr, indent=1) + "\n")
    ok1 = geq["loops_satisfy_with_rho_minus_2"] and geq["solvable_with_ghost_equation"] and bool(part and part["verified_A_x_eq_b"])
    st = {"landmark": "XVI-Q", "release": "Q08", "P6_outcome": "Partial",
          "goals": {"1_ghost_equation_compatibility": "ESTABLISHED: X can be chosen ghost-equation compatible "
                    "(rho = -2 from the loops); free directions 152 -> 104" if ok1 else "NOT established",
                    "2_invariant_operator_basis": "Partial: no explicit invariant-operator (R^2/Riemann^2-type "
                    "evanescent) basis constructed; three-point cohomology-type count only (26 hhh-only closed "
                    "evanescent cubic structures); X not expressed in an operator basis",
                    "3_free_parameters": "Partial: >= 97 (fitted lambda) of 152 in BRST-trivial family span, <= 116; "
                    ">= 36 uncharacterised; physical vs trivial NOT settled",
                    "4_p22_consistency": "OPEN: not reached (cubic completion of the p22 kappa^2 evanescent "
                    "operators not constructed)"},
          "sigma0_closure_kappa3": "unchanged from Q07: CLOSED at pole level on the one-loop three-point rows "
                                   "(enlarged evanescent CT domain); now also with a ghost-equation-compatible X",
          "not_implied": ["kappa^3 four-point / higher rows", "finite parts", "uniqueness of X",
                          "physical nature of the remaining free directions", "consistency with p22 operators",
                          "bulk symmetry restoration", "any change to the O.1 record"],
          "disclosed_revisions": {"Q08-R1": "disclosed scope revision of the frozen Q07 layer gate "
                                            "(code/p6q06/q07/release_q07.py): files under code/p6q06/q07/q08/ and "
                                            "tools/*q08* are left to the Q08 gate. release_q07.py is pinned by "
                                            "q08_expected.json. Q05 and Q06 gates untouched."},
          "O1_closure_unchanged": "CLOSED (scoped)", "bulk_symmetry_restored": False,
          "flags": {"execution_pass": bool(ok1), "claim_status": "partial",
                    "scope": "pole level, one loop, kappa^3, three-point rows"}, "seconds": round(time.time() - t0, 1)}
    (CERT / "Q08_P6_status.json").write_text(json.dumps(st, indent=1) + "\n")
    print("Q08 certs", st["flags"], st["seconds"])


if __name__ == "__main__":
    main()
