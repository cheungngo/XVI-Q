#!/usr/bin/env python3
"""XVI-Q Q09 certificates: Q07-review fixes (sigma0 relabel, Q04 independent report bound by hash, scan extension)
and goal (a) derived normalisations of the BRST-trivial families; (b)-(d) status."""
from __future__ import annotations
import hashlib, json, subprocess, sys, time  # noqa: E401
from pathlib import Path

HERE = Path(__file__).resolve().parent
Q = HERE.parents[4]
CERT = Q / "certificates"
LOG = Q / "build" / "logs"
sys.path.insert(0, str(HERE))
Q04_REPORT = "reproducibility/independent_rebuild_report_q04.json"
Q04_REPORT_SHA256 = "bc2b16d99f01cf172a5c67f6ee1abb5fdb7b9b5c30e6ad3e510c241ffa442a45"


def main():
    t0 = time.time()
    r = subprocess.run([sys.executable, "code/p6q06/q07/q08/q09/fam9.py"], cwd=Q, capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(r.stderr[-2000:])
    fam = json.loads((LOG / "q09_fam.json").read_text())
    import guard_scan3
    sc = guard_scan3.report()
    (CERT / "Q09_guard_scan3.json").write_text(json.dumps(sc, indent=1) + "\n")
    q4 = hashlib.sha256((Q / Q04_REPORT).read_bytes()).hexdigest()
    (CERT / "Q09_q04_report_binding.json").write_text(json.dumps(
        {"file": Q04_REPORT, "sha256": q4, "pinned": Q04_REPORT_SHA256, "bound": q4 == Q04_REPORT_SHA256}, indent=1) + "\n")
    tr = json.loads((CERT / "Q08_free_directions.json").read_text())
    fa = {"derivation": ("sigma(Y) = (S0, Y). Y = K.F: X_hhh = D-op(F), X_Khc = A-op(F); the row operators are the same "
                         "A-op and D-op and commute, so sigma^2 = 0 fixes lambda_F = -1 (relative sign of the field / "
                         "antifield derivatives of a fermionic Y); likewise lambda_P = -1 with the E-op; lambda_G = +1 "
                         "(both pieces are dY/dc* . dS0/dc)."),
          "families": {k: v for k, v in fam.items() if k in ("F", "G", "P")},
          "dim_trivial_span_derived": fam.get("dim_trivial_span_derived"),
          "dim_trivial_span_derived_inside_ghost_eq": fam.get("dim_trivial_span_derived_inside_ghost_eq"),
          "Q08_fitted_counts_superseded": {"fitted_lower_bound_in_N": tr["dim_N_cap_all_families_best_lambda"],
                                           "status": "fitted, not derivation-backed; superseded by the derived counts"},
          "diagnosis": ("F at the derived lambda: 139/204 generators in N; every one of the 65 failures carries a "
                        "trace monomial (K or h trace), pointing to the eta-trace / d-convention inside the probe "
                        "contraction rather than to the derivation. The Q08 G and P constructions are NOT closed at the "
                        "derived lambda (0/33, 0/134): their pieces were built by hand and are inconsistent with the row "
                        "operators. Both are OPEN.")}
    (CERT / "Q09_trivial_families.json").write_text(json.dumps(fa, indent=1) + "\n")
    n_g = tr["homogeneous_space"]["dim_N_with_ghost_eq"]
    st = {"landmark": "XVI-Q", "release": "Q09", "P6_outcome": "Partial",
          "sigma0_status_kappa3": "evanescent kappa^3 three-point defect sigma0-exact at pole level; X non-unique "
                                  "(ghost-eq compatible)",
          "relabel_note": "supersedes the Q07/Q08 wording 'sigma0 closure CLOSED at pole level' (Q07 review, medium); "
                          "frozen Q07/Q08 certificates left unchanged",
          "goals": {"a_derived_normalisations": "Partial: lambda_F = lambda_P = -1, lambda_G = +1 derived; F verified on "
                    "139/204 generators (failures all trace-type, OPEN); G, P constructions not closed (OPEN)",
                    "b_classification": f"Partial: derived trivial span {fam.get('dim_trivial_span_derived')} "
                    f"({fam.get('dim_trivial_span_derived_inside_ghost_eq')} inside the ghost-equation subspace of dim "
                    f"{n_g}); the remaining directions are NOT classified (trivial vs physical)",
                    "c_invariant_operator_basis": "OPEN: not constructed; X not expressed in an operator basis",
                    "d_p22_consistency": "OPEN: not reached"},
          "review_fixes": {"sigma0_relabel": True, "Q04_report_bound": q4 == Q04_REPORT_SHA256,
                           "scan_extension": sc["extension_clean"] and sc["all_probes_flagged"] and sc["base_scan2_clean"]},
          "not_implied": ["uniqueness or physicality of X", "kappa^3 four-point / higher rows", "finite parts",
                          "p22 consistency", "bulk symmetry restoration", "any change to the O.1 record"],
          "disclosed_revisions": {"Q09-R1": "disclosed scope revision of the frozen Q08 layer gate "
                                            "(code/p6q06/q07/q08/release_q08.py): files under code/p6q06/q07/q08/q09/ "
                                            "and tools/*q09* are left to the Q09 gate; release_q08.py is pinned by "
                                            "q09_expected.json. Q05-Q07 gates untouched. guard_scan2 (Q05) untouched; "
                                            "the extension is a new Q09 module."},
          "O1_closure_unchanged": "CLOSED (scoped)", "bulk_symmetry_restored": False,
          "flags": {"execution_pass": True, "claim_status": "partial",
                    "scope": "pole level, one loop, kappa^3, three-point rows"}, "seconds": round(time.time() - t0, 1)}
    (CERT / "Q09_P6_status.json").write_text(json.dumps(st, indent=1) + "\n")
    print("Q09 certs", st["review_fixes"], st["seconds"])


if __name__ == "__main__":
    main()
