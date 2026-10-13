#!/usr/bin/env python3
"""XVI-Q Q02: runs the bulk graphs (bulk_graphs.py) and the bulk identity / O.1 walk (bulk_identity.py) and writes
certificates Q02_bulk_graphs.json, Q02_D12_reconciliation.json, Q02_bulk_identity.json, Q02_O1.json."""
from __future__ import annotations
import json
import subprocess
import sys
from pathlib import Path

Q = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Q / "code" / "q02"))
import bulk_identity as BI  # noqa: E402

CERT = Q / "certificates"


def w(name, obj):
    (CERT / f"{name}.json").write_text(json.dumps({"landmark": "XVI-Q", "release": "Q02", **obj}, indent=1, default=str) + "\n")


def main():
    bg_path = Q / "build" / "q02_bulk_graphs.json"
    r = subprocess.run([sys.executable, str(Q / "code" / "q02" / "bulk_graphs.py"), str(bg_path)], capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-2000:], r.stderr[-3000:])
        raise SystemExit("bulk_graphs failed")
    bg = json.loads(bg_path.read_text())
    w("Q02_bulk_graphs", bg)
    A = json.loads((Q / "provenance/XVI-P/certificates/P15_batchA.json").read_text())
    s23 = (Q / "provenance/XVI-P/certificates/P23_O1_scope.json").read_text()
    rec = {"graph": "D12_hh_scalar_tadpole",
           "XVI_O": "inventory 'provisional_listed', sym_factor 1/2 'Unverified', never evaluated; XVI-O ledgers do NOT label "
                    "D12 'Missing' (searched: O_quantum_target, amplitude_status, symmetry_defect, diagram_inventory)",
           "XVI_P_p14r": A["D12_slot"]["previous_status"], "XVI_P_p15": A["D12_slot"]["status"],
           "XVI_P_p15_derivation": A["D12_slot"]["derivation"],
           "XVI_P_p23_stale_line_present": "D12 Missing never zero per P14" in s23,
           "XVI_Q_Q02_independent": bg["tadpoles"]["D12"],
           "reconciliation": ("'Missing' was the p14r status (integrand only, no value) and survives only as stale text in XVI-P "
                              "P23_O1_scope's O.5 evidence line (self-contradictory with P15, already flagged by the p23 review). "
                              "D12 was EVALUATED in XVI-P p15: vanishing_under_declared_prescription by two routes with IR-mass "
                              "and cutoff analysis; Q02 confirms independently on the engine integrand: every member "
                              "int (l^2)^j / l^2 with omega in {2, 4}, no logarithmic member, so neither a UV nor an IR pole: DR "
                              "value 0 DERIVED (not by convention; Missing is never set to zero). Recorded as erratum E6; XVI-O "
                              "and XVI-P untouched."),
           "status": "Evaluated_vanishing_under_declared_prescription",
           "consistent": A["D12_slot"]["status"] == "vanishing_under_declared_prescription"
           and not bg["tadpoles"]["D12"]["logarithmic_member"]}
    w("Q02_D12_reconciliation", rec)
    sys.path.insert(0, str(Q / "code" / "q02"))
    import d04_route4 as R4
    scr = R4.main()
    scr["engine_routes"] = {k: bg["D04"][k] for k in ("route1_A", "route2_IBP_A", "routes_agree", "N_factorised",
                                                       "N_proportional_to_p2_times_scalar_denominator")}
    scr["engine_mutation_drop_trace"] = bg["mutations"].get("D04_drop_trace_part_value")
    scr["E7"] = ("symmetry factor of D04 is 1 (distinct-species internal lines; no automorphism with labelled external legs). "
                 "XVI-O's '1/2' was Unverified and never used for a value; XVI-P never evaluated D04 (excluded from ch by "
                 "phi-parity) and its D11 (two identical phi lines) keeps S_f = 1/2 correctly: NO earlier value changes.")
    scr["D05a_D12"] = ("genuine one-propagator tadpoles: no external momentum flows through the loop line, so the integral "
                       "is scaleless; omega = 2 (D05a) / 2, 4 (D12) power divergences, no logarithmic member -> no 1/eps on "
                       "either side; 0 is not a UV/IR pole cancellation. (Unlike D04, whose bubble is NOT scaleless.)")
    w("Q02_D04_scrutiny", scr)
    bulk = BI.main(bg_path)
    w("Q02_bulk_identity", bulk)
    o2 = json.loads((CERT / "Q01_O2.json").read_text())["O2_kappa2_master_status"]
    walk = BI.o1_walk(bulk, o2, rec["consistent"], scr)
    walk["scope"] = ("one-loop kappa^2 bulk two-point sector (hh, phi phi, cbar c, K'.c, one-point) of S_inv, m = 0, frozen gauge "
                     "fermion, under the XVI-P declared prescription (DRED with opaque d, MS-bar, Feynman +i0, p22 evanescent "
                     "counterterms; external momenta in QdS); insertion frozen (XVI-N massless obstruction untouched)")
    walk["not_implied"] = ["bulk_symmetry_restored (all orders / all n-point): false",
                           "finite renormalised amplitudes with non-QdS (epsilon-scalar) external components: the p22 "
                           "Q4S-trace O(eps)x pole remnant stays disclosed (identity unaffected: projected Q4S sum 0)",
                           "three-point / higher orders (P.6), prepared S_D, insertion"]
    w("Q02_O1", walk)
    print(json.dumps({"d04_declared_zero": scr["declared_zero"], "bulk": bulk["bulk_identity_outcome"], "rows": bulk["n_rows"], "unclassified": bulk["unclassified_rows"],
                      "O1": walk["O1_closure"], "remaining": walk["remaining"]}))


if __name__ == "__main__":
    main()
