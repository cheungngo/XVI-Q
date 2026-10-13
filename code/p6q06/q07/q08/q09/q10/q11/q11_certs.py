#!/usr/bin/env python3
"""XVI-Q Q11 certificates: sigma-images in the gauge-fixed BV variables (Ktilde = K + cbar.F'), diagnosis of the
Q08 F defect, conditional classification of the ghost-equation free space N_g, status."""
from __future__ import annotations
import json, subprocess, sys, time  # noqa: E401
from pathlib import Path

HERE = Path(__file__).resolve().parent
Q = HERE.parents[6]
CERT = Q / "certificates"
LOG = Q / "build" / "logs"
STEPS = ["diag11.py", "sigma11.py", "scan11.py", "trace11.py", "class11.py"]


def main():
    t0 = time.time()
    for s in STEPS:
        r = subprocess.run([sys.executable, str(HERE / s)], cwd=Q, capture_output=True, text=True)
        if r.returncode:
            raise RuntimeError(s + r.stderr[-2000:])
    diag = json.loads((LOG / "q11_diag.json").read_text())
    sig = json.loads((LOG / "q11_sigma.json").read_text())
    scan = json.loads((LOG / "q11_scan.json").read_text())
    tr = json.loads((LOG / "q11_trace.json").read_text())
    cl = json.loads((LOG / "q11_class.json").read_text())
    sig.pop("seconds", None)
    si = {"construction": ("Y built from Ktilde = K + cbar.F'(h) (canonical shift by Psi = cbar.F, B integrated out); "
                           "the cbar partner of every K-functional is rho*T(.), rho = -2 (Q08 ghost equation). "
                           "sigma = (S0, .) applied with the row channel operators (D, A, E), antibracket sign -1 "
                           "between antifield-type and field-type pieces; G family includes the c*cc piece "
                           "(h -> s0 h, antisymmetrised in the two c's)."),
          "Q08_F_defect_diagnosis": {"blocks_violated_by_Q08_pieces": diag,
                                     "reading": ("the Q08 F pieces violate only the h.h.c row; 12 generators fail through "
                                                 "the antifield-type piece alone, i.e. Q08 F omitted the cbar partner "
                                                 "(Ktilde) and the E channel")},
          "results_at_derived_signs": sig,
          "lam_mu_scan_diagnostic": scan,
          "trace_correlation": tr,
          "reading": ("With the Ktilde construction 71 of 194 nonzero F images and 8 of 33 G images lie in N (all of them "
                      "also satisfy the ghost equation). The remaining failures are all in the h.h.c row (G also in the "
                      "ghost rows). No single normalisation works for every generator (lam/mu scan), and failures occur "
                      "with and without trace monomials, so the eta-trace hypothesis is NOT the (sole) cause. The "
                      "row operators as inserted on evanescent X do not compose to zero with this sigma: sigma^2 = 0 "
                      "is not realised on the evanescent sector of the h.h.c row. Diagnosis: OPEN (candidate: the "
                      "4D-momentum rule T(p, .) = 0 for hat components, under which s0 h_hat = 0, versus the "
                      "d-dimensional gauge structure of the kinetic operator).")}
    (CERT / "Q11_sigma_images.json").write_text(json.dumps(si, indent=1) + "\n")
    cla = {"N_g": cl["dim_N_g"], "sigma_span_passing": cl["dim_sigma_span_passing"],
           "blind_in_N_g": cl["dim_blind_in_N_g"], "sigma_cap_blind": cl["dim_sigma_cap_blindg"],
           "sigma_plus_blind": cl["dim_sigma_plus_blindg"], "unclassified": cl["dim_remaining_N_g_mod_(sigma+blind)"],
           "status": ("CONDITIONAL: the 46 directions are sigma-images of local Ktilde-functionals that lie in N_g, hence "
                      "trivial IF this sigma is the correct linearised Slavnov-Taylor operator on the evanescent sector; "
                      "that is exactly what the failures above put in doubt. 24 blind directions (9 overlapping the "
                      "sigma span). 43 directions are unclassified: NOT shown to be physical."),
           "physical_evanescent_directions": "not established (no direction certified physical)"}
    (CERT / "Q11_classification.json").write_text(json.dumps(cla, indent=1) + "\n")
    st = {"landmark": "XVI-Q", "release": "Q11", "P6_outcome": "Partial",
          "sigma0_status_kappa3": "evanescent kappa^3 three-point defect sigma0-exact at pole level; X non-unique "
                                  "(ghost-eq compatible)",
          "goals": {"1_sigma_images": "Partial: genuine Ktilde construction built (F, G; P absorbed as the cbar partner of "
                                      "F); 71/194 F and 8/33 G in N; Q08 F defect diagnosed (missing cbar partner / E "
                                      "channel); remaining failures OPEN (not trace-driven)",
                    "2_classification": "Partial (conditional): N_g 104 = sigma span 46 + blind 24 - overlap 9 + 43 "
                                        "unclassified",
                    "3_invariant_basis": "OPEN: not constructed",
                    "3_p22_consistency": "OPEN: not reached"},
          "not_implied": ["uniqueness or physicality of X", "an unconditional count of BRST-trivial directions",
                          "any physical evanescent direction", "p22 consistency", "bulk symmetry restoration",
                          "any change to the O.1 record"],
          "disclosed_revisions": {"Q11-R1": "disclosed scope revision of the frozen Q10 layer gate "
                                            "(code/p6q06/q07/q08/q09/q10/release_q10.py): files under "
                                            "code/p6q06/q07/q08/q09/q10/q11/ and tools/*q11* are left to the Q11 gate; "
                                            "and the Q10 gate treats n_files_scanned of Q09_guard_scan3 as a lower bound (new Q11 files are scanned, not exempted); release_q10.py is pinned by q11_expected.json. Q05-Q09 gates untouched."},
          "premises_note": "OpenAI catalogue and En05-En07 sheets are not premises (Q10 intake)",
          "O1_closure_unchanged": "CLOSED (scoped)", "bulk_symmetry_restored": False, "insertion": "frozen",
          "flags": {"execution_pass": True, "claim_status": "partial",
                    "scope": "pole level, one loop, kappa^3, three-point rows"}, "seconds": round(time.time() - t0, 1)}
    (CERT / "Q11_P6_status.json").write_text(json.dumps(st, indent=1) + "\n")
    print("Q11 certs", st["seconds"])


if __name__ == "__main__":
    main()
